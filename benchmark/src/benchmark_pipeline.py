from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

try:
    import igl
    import igl.copyleft.cgal as igl_cgal
except ImportError as exc:
    raise RuntimeError(
        "Python package 'libigl' is required. Install dependencies from benchmark/requirements.txt."
    ) from exc

try:
    import psutil
except ImportError:
    psutil = None


SUPPORTED_MESH_SUFFIXES = {".obj", ".off", ".stl", ".ply"}
MEMORY_CHECK_INTERVAL_SEC = 0.01
OP_MAP = {
    "union": "union",
    "u": "union",
    "intersection": "intersection",
    "intersect": "intersection",
    "i": "intersection",
    "difference": "difference",
    "minus": "difference",
    "m": "difference",
    "d": "difference",
}
OP_TOKEN = {"union": "U", "intersection": "I", "difference": "D"}


def _normalize_op(operation: str) -> str:
    key = operation.strip().lower()
    if key not in OP_MAP:
        raise ValueError(f"Unsupported operation: {operation}")
    return OP_MAP[key]


def _to_serializable(obj: Any) -> Any:
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(v) for v in obj]
    return obj


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_to_serializable(payload), handle, indent=2)


def _write_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with path.open("w", encoding="utf-8", newline="") as handle:
            handle.write("")
        return
    headers = sorted({k for row in rows for k in row.keys()})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _to_serializable(row.get(k)) for k in headers})


def load_mesh(mesh_path: Path) -> Tuple[np.ndarray, np.ndarray]:
    v, f = igl.read_triangle_mesh(str(mesh_path))
    v = np.asarray(v, dtype=np.float64)
    f = np.asarray(f, dtype=np.int64)
    if v.size == 0 and f.size == 0:
        return np.zeros((0, 3), dtype=np.float64), np.zeros((0, 3), dtype=np.int64)
    if v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or f.shape[1] != 3:
        raise ValueError(f"Expected triangle mesh with 3D vertices: {mesh_path}")
    return v, f


def save_mesh(mesh_path: Path, v: np.ndarray, f: np.ndarray) -> None:
    mesh_path.parent.mkdir(parents=True, exist_ok=True)
    v = np.asarray(v, dtype=np.float64)
    f = np.asarray(f, dtype=np.int64)
    if v.size == 0 or f.size == 0:
        with mesh_path.open("w", encoding="utf-8") as handle:
            handle.write("# empty mesh\n")
        return
    igl.write_triangle_mesh(str(mesh_path), v, f)


def remove_unreferenced(v: np.ndarray, f: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if f.size == 0:
        return np.zeros((0, 3), dtype=np.float64), np.zeros((0, 3), dtype=np.int64)
    used = np.unique(f.reshape(-1))
    remap = np.full(v.shape[0], -1, dtype=np.int64)
    remap[used] = np.arange(used.shape[0], dtype=np.int64)
    v_new = v[used]
    f_new = remap[f]
    return v_new, f_new


def signed_volume(v: np.ndarray, f: np.ndarray) -> float:
    if f.size == 0:
        return 0.0
    tri = v[f]
    return float(np.sum(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))) / 6.0)


def face_areas(v: np.ndarray, f: np.ndarray) -> np.ndarray:
    if f.size == 0:
        return np.zeros((0,), dtype=np.float64)
    tri = v[f]
    return 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)


def boundary_edges(f: np.ndarray) -> np.ndarray:
    if f.size == 0:
        return np.zeros((0, 2), dtype=np.int64)
    e = np.vstack((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]])).astype(np.int64)
    e_sorted = np.sort(e, axis=1)
    uniq, counts = np.unique(e_sorted, axis=0, return_counts=True)
    return uniq[counts == 1]


def mesh_bbox(v: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if v.size == 0:
        return np.zeros(3), np.zeros(3)
    return v.min(axis=0), v.max(axis=0)


def mesh_stats(v: np.ndarray, f: np.ndarray) -> Dict[str, Any]:
    bmin, bmax = mesh_bbox(v)
    return {
        "num_vertices": int(v.shape[0]),
        "num_faces": int(f.shape[0]),
        "bbox_min": bmin.tolist(),
        "bbox_max": bmax.tolist(),
        "bbox_diag": float(np.linalg.norm(bmax - bmin)),
        "boundary_edges": int(boundary_edges(f).shape[0]),
        "surface_area": float(face_areas(v, f).sum()),
    }


def uv_sphere_mesh(
    center: np.ndarray,
    radius: float,
    stacks: int,
    slices: int,
) -> Tuple[np.ndarray, np.ndarray]:
    stacks = max(3, int(stacks))
    slices = max(6, int(slices))

    verts: List[List[float]] = [[0.0, 0.0, 1.0]]
    for i in range(1, stacks):
        phi = math.pi * i / stacks
        sp = math.sin(phi)
        cp = math.cos(phi)
        for j in range(slices):
            theta = 2.0 * math.pi * j / slices
            verts.append([sp * math.cos(theta), sp * math.sin(theta), cp])
    verts.append([0.0, 0.0, -1.0])
    bottom = len(verts) - 1

    def ring(i: int, j: int) -> int:
        return 1 + (i - 1) * slices + (j % slices)

    faces: List[List[int]] = []
    for j in range(slices):
        faces.append([0, ring(1, j), ring(1, j + 1)])
    for i in range(1, stacks - 1):
        for j in range(slices):
            a = ring(i, j)
            b = ring(i, j + 1)
            c = ring(i + 1, j)
            d = ring(i + 1, j + 1)
            faces.append([a, c, b])
            faces.append([b, c, d])
    for j in range(slices):
        faces.append([ring(stacks - 1, j + 1), ring(stacks - 1, j), bottom])

    v = np.asarray(verts, dtype=np.float64)
    f = np.asarray(faces, dtype=np.int64)
    if signed_volume(v, f) < 0:
        f = f[:, [0, 2, 1]]
    v = v * float(radius) + center.reshape(1, 3)
    return v, f


def index_dataset(dataset_dir: Path) -> List[Dict[str, Any]]:
    mesh_paths = sorted(
        p for p in dataset_dir.rglob("*") if p.is_file() and p.suffix.lower() in SUPPORTED_MESH_SUFFIXES
    )
    records: List[Dict[str, Any]] = []
    for idx, mesh_path in enumerate(mesh_paths):
        record: Dict[str, Any] = {
            "mesh_id": f"mesh_{idx:05d}",
            "path": str(mesh_path),
            "suffix": mesh_path.suffix.lower(),
            "status": "ok",
        }
        try:
            v, f = load_mesh(mesh_path)
            record.update(mesh_stats(v, f))
        except Exception as exc:  # explicit capture for manifesting failures
            record["status"] = "error"
            record["error"] = str(exc)
        records.append(record)
    return records


def _edge_manifold_bool(f: np.ndarray) -> bool:
    edge_raw = igl.is_edge_manifold(f)
    if isinstance(edge_raw, tuple):
        return bool(edge_raw[0])
    if isinstance(edge_raw, np.ndarray):
        return bool(np.all(edge_raw))
    return bool(edge_raw)


def _vertex_manifold_bool(f: np.ndarray) -> bool:
    vertex_raw = igl.is_vertex_manifold(f)
    if isinstance(vertex_raw, tuple):
        return bool(np.all(np.asarray(vertex_raw[0])))
    if isinstance(vertex_raw, np.ndarray):
        return bool(np.all(vertex_raw))
    return bool(vertex_raw)


def _topology_faces_for_checks(v: np.ndarray, f: np.ndarray, suffix: str) -> Tuple[np.ndarray, bool]:
    if f.size == 0:
        return f, False
    if suffix != ".stl":
        return f, False
    bmin, bmax = mesh_bbox(v)
    diag = float(np.linalg.norm(bmax - bmin))
    epsilon = max(diag * 1e-12, 1e-15)
    sv, _, _, sf = igl.remove_duplicate_vertices(v, f, epsilon)
    sf = np.asarray(sf, dtype=np.int64)
    if sv.size == 0 or sf.size == 0:
        return f, False
    return sf, True


def compute_dataset_mesh_stats(dataset_dir: Path) -> List[Dict[str, Any]]:
    mesh_paths = sorted(
        p for p in dataset_dir.rglob("*") if p.is_file() and p.suffix.lower() in SUPPORTED_MESH_SUFFIXES
    )
    records: List[Dict[str, Any]] = []
    for idx, mesh_path in enumerate(mesh_paths):
        size_bytes = int(mesh_path.stat().st_size)
        record: Dict[str, Any] = {
            "mesh_id": f"mesh_{idx:05d}",
            "path": str(mesh_path),
            "suffix": mesh_path.suffix.lower(),
            "size_bytes": size_bytes,
            "size_mb": float(size_bytes / (1024.0 * 1024.0)),
            "status": "ok",
            "num_triangles": 0,
            "is_closed": False,
            "is_manifold": False,
            "is_edge_manifold": False,
            "is_vertex_manifold": False,
            "topology_faces_source": "raw",
        }
        try:
            v, f = load_mesh(mesh_path)
            topo_f, welded_for_topology = _topology_faces_for_checks(v, f, mesh_path.suffix.lower())
            boundary_count = int(boundary_edges(topo_f).shape[0])
            is_closed = boundary_count == 0
            is_edge_manifold = _edge_manifold_bool(topo_f)
            is_vertex_manifold = _vertex_manifold_bool(topo_f)
            is_manifold = bool(is_edge_manifold and is_vertex_manifold)
            record.update(mesh_stats(v, f))
            record["num_triangles"] = int(f.shape[0])
            record["is_closed"] = bool(is_closed)
            record["is_manifold"] = bool(is_manifold)
            record["is_edge_manifold"] = bool(is_edge_manifold)
            record["is_vertex_manifold"] = bool(is_vertex_manifold)
            record["topology_faces_source"] = "welded" if welded_for_topology else "raw"
        except Exception as exc:
            record["status"] = "error"
            record["error"] = str(exc)
        records.append(record)
    return records


def get_dataset_mesh_stats(dataset_dir: Path, update_stats: bool = False) -> List[Dict[str, Any]]:
    stats_json = dataset_dir / "mesh_stats.json"
    stats_csv = dataset_dir / "mesh_stats.csv"
    if stats_json.exists() and not update_stats:
        with stats_json.open("r", encoding="utf-8") as handle:
            loaded = json.load(handle)
        if isinstance(loaded, list):
            has_required_fields = all(
                isinstance(rec, dict)
                and "is_closed" in rec
                and "is_manifold" in rec
                and "num_triangles" in rec
                and "topology_faces_source" in rec
                for rec in loaded
            )
            if has_required_fields:
                return loaded
    records = compute_dataset_mesh_stats(dataset_dir)
    _write_json(stats_json, records)
    _write_csv(stats_csv, records)
    return records


def sample_random_cases(
    mesh_records: Sequence[Dict[str, Any]],
    num_pairs: int,
    operations: Sequence[str],
    seed: int,
) -> List[Dict[str, Any]]:
    valid = [r for r in mesh_records if r.get("status") == "ok"]
    if not operations:
        raise ValueError("At least one operation must be provided.")
    if len(valid) < 2:
        return []

    rng = np.random.default_rng(seed)
    normalized_ops = [_normalize_op(op) for op in operations]
    unique_pairs = [(i, j) for i in range(len(valid)) for j in range(i + 1, len(valid))]
    if not unique_pairs:
        return []
    rng.shuffle(unique_pairs)
    pair_count = min(int(num_pairs), len(unique_pairs))
    cases: List[Dict[str, Any]] = []
    for i, (ia, ib) in enumerate(unique_pairs[:pair_count]):
        op = normalized_ops[int(rng.integers(0, len(normalized_ops)))]
        cases.append(
            {
                "case_id": f"case_{i:05d}",
                "input_a": valid[ia]["path"],
                "input_b": valid[ib]["path"],
                "operation": op,
                "seed_spheres_a": int(rng.integers(0, 2**31 - 1)),
                "seed_spheres_b": int(rng.integers(0, 2**31 - 1)),
            }
        )
    return cases


def sample_cut_spheres(
    v: np.ndarray,
    f: np.ndarray,
    sphere_cfg: Dict[str, Any],
    rng_seed: int,
) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(rng_seed)
    bmin, bmax = mesh_bbox(v)
    diag = float(np.linalg.norm(bmax - bmin))
    n_min = int(sphere_cfg["count_min"])
    n_max = int(sphere_cfg["count_max"])
    radius_min = float(sphere_cfg["radius_min_frac"]) * diag
    radius_max = float(sphere_cfg["radius_max_frac"]) * diag
    center_mode = str(sphere_cfg.get("center_mode", "surface")).lower()
    stacks = int(sphere_cfg["stacks"])
    slices = int(sphere_cfg["slices"])

    count = int(rng.integers(n_min, n_max + 1))
    spheres: List[Dict[str, Any]] = []
    for i in range(count):
        if center_mode == "bbox":
            center = rng.uniform(bmin, bmax).astype(np.float64)
        else:
            if f.shape[0] == 0:
                center = rng.uniform(bmin, bmax).astype(np.float64)
                radius = float(rng.uniform(radius_min, radius_max))
                sv, sf = uv_sphere_mesh(center, radius, stacks, slices)
                spheres.append(
                    {
                        "sphere_id": i,
                        "center": center,
                        "radius": radius,
                        "v": sv,
                        "f": sf,
                    }
                )
                continue
            fi = int(rng.integers(0, f.shape[0]))
            tri = v[f[fi]]
            r1 = float(rng.uniform(0.0, 1.0))
            r2 = float(rng.uniform(0.0, 1.0))
            s = math.sqrt(r1)
            w0 = 1.0 - s
            w1 = s * (1.0 - r2)
            w2 = s * r2
            center = (w0 * tri[0] + w1 * tri[1] + w2 * tri[2]).astype(np.float64)
        radius = float(rng.uniform(radius_min, radius_max))
        sv, sf = uv_sphere_mesh(center, radius, stacks, slices)
        spheres.append(
            {
                "sphere_id": i,
                "center": center,
                "radius": radius,
                "v": sv,
                "f": sf,
            }
        )
    return spheres


def apply_sphere_cuts(
    v: np.ndarray,
    f: np.ndarray,
    provenance: np.ndarray,
    spheres: Sequence[Dict[str, Any]],
    target_tag: Optional[int],
    save_steps_dir: Optional[Path],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    v_cur = np.asarray(v, dtype=np.float64)
    f_cur = np.asarray(f, dtype=np.int64)
    prov_cur = np.asarray(provenance, dtype=np.int64)
    cut_stats: List[Dict[str, Any]] = []

    if save_steps_dir is not None:
        save_steps_dir.mkdir(parents=True, exist_ok=True)

    for cut_idx, sphere in enumerate(spheres):
        if f_cur.size == 0:
            cut_stats.append(
                {
                    "cut_index": cut_idx,
                    "removed_faces": 0,
                    "remaining_faces": 0,
                    "removed_fraction": 0.0,
                    "skipped_empty_input": True,
                }
            )
            continue

        vt, ft, d, j = igl_cgal.trim_with_solid(v_cur, f_cur, sphere["v"], sphere["f"])
        d = np.asarray(d).reshape(-1).astype(bool)
        j = np.asarray(j).reshape(-1).astype(np.int64)

        parent_tags = prov_cur[j]
        # trim_with_solid's boolean label convention in Python bindings may vary
        # across versions; remove the minority side to realize local "hole" cuts.
        remove_base = d
        if remove_base.sum() > (remove_base.size - remove_base.sum()):
            remove_base = ~remove_base
        if target_tag is None:
            remove_mask = remove_base
        else:
            remove_mask = remove_base & (parent_tags == int(target_tag))
        keep_mask = ~remove_mask

        f_keep = np.asarray(ft, dtype=np.int64)[keep_mask]
        prov_keep = parent_tags[keep_mask]
        removed_faces = int(remove_mask.sum())
        total_faces = int(ft.shape[0])
        removed_fraction = (removed_faces / total_faces) if total_faces else 0.0

        v_cur, f_cur = remove_unreferenced(np.asarray(vt, dtype=np.float64), f_keep)
        prov_cur = prov_keep

        if save_steps_dir is not None:
            save_mesh(save_steps_dir / f"step_{cut_idx:02d}.obj", v_cur, f_cur)

        cut_stats.append(
            {
                "cut_index": cut_idx,
                "removed_faces": removed_faces,
                "remaining_faces": int(f_cur.shape[0]),
                "removed_fraction": removed_fraction,
                "skipped_empty_input": False,
            }
        )

    return v_cur, f_cur, prov_cur, cut_stats


def run_boolean(v_a: np.ndarray, f_a: np.ndarray, v_b: np.ndarray, f_b: np.ndarray, operation: str):
    op = _normalize_op(operation)
    v_x, f_x, j_x = igl_cgal.mesh_boolean(v_a, f_a, v_b, f_b, op)
    return np.asarray(v_x, dtype=np.float64), np.asarray(f_x, dtype=np.int64), np.asarray(j_x).reshape(-1).astype(np.int64)


def _operation_argument(operation: str, operation_map: Dict[str, str]) -> str:
    op = _normalize_op(operation)
    key_candidates = {
        "union": ["union", "u"],
        "intersection": ["intersection", "intersect", "int", "i"],
        "difference": ["difference", "diff", "minus", "m", "d"],
    }[op]
    lowered = {str(k).lower(): str(v) for k, v in operation_map.items()}
    for key in key_candidates:
        if key in lowered:
            return lowered[key]
    raise ValueError(f"operation_map does not define mapping for '{op}'.")


def _find_new_mesh_output(
    working_directory: Path,
    start_time: float,
    exclude_files: Sequence[Path],
) -> Optional[Path]:
    exclude_resolved = {p.resolve() for p in exclude_files if p.exists()}
    candidates: List[Tuple[float, Path]] = []
    for suffix in SUPPORTED_MESH_SUFFIXES:
        for mesh_path in working_directory.rglob(f"*{suffix}"):
            try:
                resolved = mesh_path.resolve()
                if resolved in exclude_resolved:
                    continue
                stat = mesh_path.stat()
            except OSError:
                continue
            if stat.st_mtime >= start_time:
                candidates.append((stat.st_mtime, mesh_path))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def _terminate_process_tree(pid: int) -> None:
    if psutil is None:
        return
    try:
        root_proc = psutil.Process(pid)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return
    try:
        children = root_proc.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        children = []
    for child in children:
        try:
            child.terminate()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    try:
        root_proc.terminate()
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass
    alive: List[Any] = []
    try:
        _, alive = psutil.wait_procs(children + [root_proc], timeout=2.0)
    except Exception:
        alive = []
    for proc in alive:
        try:
            proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue


def _run_subprocess_with_limits(
    argv: Sequence[str],
    cwd: Path,
    timeout_limit_sec: float,
    memory_limit_mb: float,
    check_interval_sec: float,
) -> Tuple[subprocess.CompletedProcess, Optional[float], float, Optional[str]]:
    start = time.perf_counter()
    with tempfile.TemporaryFile(mode="w+", encoding="utf-8") as stdout_file, tempfile.TemporaryFile(
        mode="w+", encoding="utf-8"
    ) as stderr_file:
        popen = subprocess.Popen(
            list(argv),
            cwd=str(cwd),
            stdout=stdout_file,
            stderr=stderr_file,
            text=True,
        )

        peak_rss_bytes = 0
        terminated_reason: Optional[str] = None
        ps_proc = psutil.Process(popen.pid) if psutil is not None else None

        if ps_proc is not None:
            try:
                peak_rss_bytes = int(ps_proc.memory_info().rss)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                peak_rss_bytes = 0

        while popen.poll() is None:
            elapsed = time.perf_counter() - start
            if timeout_limit_sec > 0 and elapsed > timeout_limit_sec:
                terminated_reason = "timeout_limit_exceeded"
                try:
                    popen.terminate()
                except OSError:
                    pass
                _terminate_process_tree(popen.pid)
                break
            if ps_proc is not None:
                try:
                    rss = int(ps_proc.memory_info().rss)
                    peak_rss_bytes = max(peak_rss_bytes, rss)
                    if memory_limit_mb > 0 and (rss / (1024.0 * 1024.0)) > memory_limit_mb:
                        terminated_reason = "memory_limit_exceeded"
                        try:
                            popen.terminate()
                        except OSError:
                            pass
                        _terminate_process_tree(popen.pid)
                        break
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            time.sleep(check_interval_sec if check_interval_sec > 0 else MEMORY_CHECK_INTERVAL_SEC)

        if popen.poll() is None:
            try:
                popen.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    popen.kill()
                except OSError:
                    pass
                try:
                    popen.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    pass

        stdout_file.flush()
        stderr_file.flush()
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read()
        stderr = stderr_file.read()

    elapsed_local = float(time.perf_counter() - start)
    peak_rss_mb = (float(peak_rss_bytes) / (1024.0 * 1024.0)) if peak_rss_bytes > 0 else None
    completed = subprocess.CompletedProcess(list(argv), popen.returncode, stdout, stderr)
    return completed, peak_rss_mb, elapsed_local, terminated_reason


def run_method_under_test(
    method_cfg: Dict[str, Any],
    operation: str,
    case_dir: Path,
    input_a_path: Path,
    input_b_path: Path,
    output_path: Path,
    limits_cfg: Optional[Dict[str, Any]] = None,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[str, Any]]:
    op = _normalize_op(operation)
    t0 = time.perf_counter()
    exe_path = method_cfg.get("exe_path")
    operation_map = method_cfg.get("operation_map")
    if not isinstance(exe_path, str) or not exe_path.strip():
        raise ValueError("Method config must provide non-empty 'exe_path'.")
    if not isinstance(operation_map, dict):
        raise ValueError("Method config must provide 'operation_map' dictionary.")

    in_a = input_a_path
    in_b = input_b_path
    out_y = output_path
    if out_y.exists():
        out_y.unlink()
    input_size_bytes = 0
    for p in (in_a, in_b):
        try:
            input_size_bytes += int(p.stat().st_size)
        except OSError:
            continue

    op_arg = _operation_argument(op, operation_map)
    base_args = [exe_path, op_arg, str(in_a), str(in_b)]
    args_with_out = [*base_args, str(out_y)]
    args_without_out = list(base_args)
    working_directory = case_dir
    exclude = [
        in_a,
        in_b,
        case_dir / "A.obj",
        case_dir / "B.obj",
        out_y,
        case_dir / "X.obj",
        case_dir / "Z.obj",
    ]
    attempts: List[Dict[str, Any]] = []
    limits_cfg = limits_cfg or {}
    timeout_limit_sec = float(limits_cfg.get("method_timeout_sec", 0.0))
    memory_limit_mb = float(limits_cfg.get("method_memory_limit_mb", 0.0))
    check_interval = float(limits_cfg.get("memory_check_interval_sec", MEMORY_CHECK_INTERVAL_SEC))

    def run_attempt(argv: List[str], label: str) -> subprocess.CompletedProcess:
        start_attempt = time.time()
        proc, peak_rss_mb, attempt_runtime, terminated_reason = _run_subprocess_with_limits(
            argv=argv,
            cwd=working_directory,
            timeout_limit_sec=timeout_limit_sec,
            memory_limit_mb=memory_limit_mb,
            check_interval_sec=check_interval,
        )
        attempts.append(
            {
                "label": label,
                "argv": argv,
                "status_code": proc.returncode,
                "stdout": proc.stdout,
                "stderr": proc.stderr,
                "runtime_sec": float(attempt_runtime),
                "peak_rss_mb": peak_rss_mb,
                "terminated_reason": terminated_reason,
            }
        )
        if not out_y.exists():
            auto_out = _find_new_mesh_output(working_directory, start_attempt, exclude)
            if auto_out is not None and auto_out.resolve() != out_y.resolve():
                shutil.copyfile(auto_out, out_y)
                attempts[-1]["auto_detected_output"] = str(auto_out)
        return proc

    proc_first = run_attempt(args_with_out, "with_output_arg")
    if (proc_first.returncode != 0 or not out_y.exists()) and not out_y.exists():
        _ = run_attempt(args_without_out, "without_output_arg")

    elapsed = time.perf_counter() - t0
    meta = {
        "runtime_sec": elapsed,
        "working_directory": str(working_directory),
        "attempts": attempts,
        "command_used": attempts[-1]["argv"] if attempts else None,
        "peak_rss_mb": max((a.get("peak_rss_mb") for a in attempts if a.get("peak_rss_mb") is not None), default=None),
        "input_size_bytes": int(input_size_bytes),
        "input_size_mb": float(input_size_bytes / (1024.0 * 1024.0)),
        "method_timeout_sec": timeout_limit_sec if timeout_limit_sec > 0 else None,
        "method_memory_limit_mb": memory_limit_mb if memory_limit_mb > 0 else None,
    }
    if not attempts:
        meta["status"] = "error"
        meta["error"] = "No process attempts executed."
        return None, None, meta

    last = attempts[-1]
    terminated_reason = last.get("terminated_reason")
    meta["status_code"] = last["status_code"]
    if last.get("terminated_reason"):
        meta["terminated_reason"] = last.get("terminated_reason")
    if meta["status_code"] != 0 or not out_y.exists():
        meta["status"] = "error"
        if last.get("terminated_reason"):
            meta["error"] = str(last.get("terminated_reason"))
        return None, None, meta
    v_y, f_y = load_mesh(out_y)
    if not bool(method_cfg.get("allow_empty_output_mesh", False)) and (v_y.shape[0] == 0 or f_y.shape[0] == 0):
        meta["status"] = "error"
        meta["error"] = "Method produced an empty output mesh."
        return None, None, meta
    meta["status"] = "ok"
    return v_y, f_y, meta


def sample_points_on_mesh(v: np.ndarray, f: np.ndarray, sample_count: int) -> np.ndarray:
    if v.size == 0 or f.size == 0 or sample_count <= 0:
        return np.zeros((0, 3), dtype=np.float64)
    _, _, p = igl.random_points_on_mesh(int(sample_count), v, f)
    return np.asarray(p, dtype=np.float64)


def _pointset_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    if a.size == 0 or b.size == 0:
        return np.zeros((0,), dtype=np.float64)
    tree = cKDTree(b)
    d, _ = tree.query(a, k=1)
    return np.asarray(d, dtype=np.float64)


def boundary_fscore(
    v_y: np.ndarray,
    f_y: np.ndarray,
    v_z: np.ndarray,
    f_z: np.ndarray,
    tol: float,
) -> Dict[str, Any]:
    e_y = boundary_edges(f_y)
    e_z = boundary_edges(f_z)

    if e_y.shape[0] == 0 and e_z.shape[0] == 0:
        return {"precision": 1.0, "recall": 1.0, "fscore": 1.0}
    if e_y.shape[0] == 0 or e_z.shape[0] == 0:
        return {"precision": 0.0, "recall": 0.0, "fscore": 0.0}

    p_y = v_y[np.unique(e_y.reshape(-1))]
    p_z = v_z[np.unique(e_z.reshape(-1))]
    d_yz = _pointset_distances(p_y, p_z)
    d_zy = _pointset_distances(p_z, p_y)
    precision = float(np.mean(d_yz <= tol))
    recall = float(np.mean(d_zy <= tol))
    fscore = float((2.0 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0)
    return {"precision": precision, "recall": recall, "fscore": fscore}


def evaluate_metrics(
    v_y: np.ndarray,
    f_y: np.ndarray,
    v_z: np.ndarray,
    f_z: np.ndarray,
    metrics_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    samples = int(metrics_cfg.get("sample_count", 8000))
    tol_frac = float(metrics_cfg.get("boundary_fscore_tol_frac", 0.001))

    py = sample_points_on_mesh(v_y, f_y, samples)
    pz = sample_points_on_mesh(v_z, f_z, samples)
    d_yz = _pointset_distances(py, pz)
    d_zy = _pointset_distances(pz, py)

    if d_yz.size == 0 or d_zy.size == 0:
        geo = {
            "hausdorff": float("nan"),
            "chamfer": float("nan"),
            "d95_y_to_z": float("nan"),
            "d95_z_to_y": float("nan"),
        }
    else:
        geo = {
            "hausdorff": float(max(d_yz.max(initial=0.0), d_zy.max(initial=0.0))),
            "chamfer": float(d_yz.mean() + d_zy.mean()),
            "d95_y_to_z": float(np.percentile(d_yz, 95)),
            "d95_z_to_y": float(np.percentile(d_zy, 95)),
        }

    bmin, bmax = mesh_bbox(v_z if v_z.size else v_y)
    tol = tol_frac * float(np.linalg.norm(bmax - bmin))
    bnd = boundary_fscore(v_y, f_y, v_z, f_z, tol)

    return {
        "geometry": geo,
        "boundary": bnd,
        "mesh_y": mesh_stats(v_y, f_y),
        "mesh_z": mesh_stats(v_z, f_z),
        "boundary_tolerance": tol,
    }


def prepare_case_data(case: Dict[str, Any], config: Dict[str, Any], prep_dir: Path) -> Dict[str, Any]:
    case_id = str(case["case_id"])
    prepared: Dict[str, Any] = {"case_id": case_id, "status": "ok"}
    try:
        prep_dir.mkdir(parents=True, exist_ok=True)
        v_a, f_a = load_mesh(Path(case["input_a"]))
        v_b, f_b = load_mesh(Path(case["input_b"]))
        operation = _normalize_op(str(case["operation"]))

        debug_cfg = config.get("debug", {})
        save_cut_meshes = bool(debug_cfg.get("save_cut_meshes", False))
        save_expected_result = bool(debug_cfg.get("save_expected_result", False))
        debug_requested = bool(save_cut_meshes or save_expected_result)
        expected_results_dir_raw = debug_cfg.get("expected_results_dir")
        expected_case_dir = Path(str(expected_results_dir_raw)) / case_id if expected_results_dir_raw else None

        if debug_requested and expected_case_dir is not None:
            expected_case_dir.mkdir(parents=True, exist_ok=True)
            if save_cut_meshes:
                save_mesh(expected_case_dir / "A.obj", v_a, f_a)
                save_mesh(expected_case_dir / "B.obj", v_b, f_b)

        v_x, f_x, j_x = run_boolean(v_a, f_a, v_b, f_b, operation)
        source_tag = (j_x >= f_a.shape[0]).astype(np.int64)  # 0=A, 1=B

        sphere_cfg = config["spheres"]
        spheres_a = sample_cut_spheres(v_a, f_a, sphere_cfg, int(case["seed_spheres_a"]))
        spheres_b = sample_cut_spheres(v_b, f_b, sphere_cfg, int(case["seed_spheres_b"]))

        if save_cut_meshes and expected_case_dir is not None:
            for s in spheres_a:
                save_mesh(expected_case_dir / "cutters" / "A" / f"sphere_{s['sphere_id']:02d}.obj", s["v"], s["f"])
            for s in spheres_b:
                save_mesh(expected_case_dir / "cutters" / "B" / f"sphere_{s['sphere_id']:02d}.obj", s["v"], s["f"])

        prov_a = np.arange(f_a.shape[0], dtype=np.int64)
        v_c, f_c, _, cuts_a = apply_sphere_cuts(
            v_a,
            f_a,
            prov_a,
            spheres_a,
            target_tag=None,
            save_steps_dir=(expected_case_dir / "C_cut_steps") if (save_cut_meshes and expected_case_dir is not None) else None,
        )

        prov_b = np.arange(f_b.shape[0], dtype=np.int64)
        v_d, f_d, _, cuts_b = apply_sphere_cuts(
            v_b,
            f_b,
            prov_b,
            spheres_b,
            target_tag=None,
            save_steps_dir=(expected_case_dir / "D_cut_steps") if (save_cut_meshes and expected_case_dir is not None) else None,
        )

        v_z_a, f_z_a, tag_after_a, z_a_stats = apply_sphere_cuts(
            v_x,
            f_x,
            source_tag,
            spheres_a,
            target_tag=0,
            save_steps_dir=(expected_case_dir / "Z_cut_A_steps") if (save_cut_meshes and expected_case_dir is not None) else None,
        )
        v_z, f_z, tag_final, z_b_stats = apply_sphere_cuts(
            v_z_a,
            f_z_a,
            tag_after_a,
            spheres_b,
            target_tag=1,
            save_steps_dir=(expected_case_dir / "Z_cut_B_steps") if (save_cut_meshes and expected_case_dir is not None) else None,
        )
        _ = tag_final

        if save_cut_meshes and expected_case_dir is not None:
            save_mesh(expected_case_dir / "X.obj", v_x, f_x)
            np.save(expected_case_dir / "X_birth_indices.npy", j_x)
            save_mesh(expected_case_dir / "C.obj", v_c, f_c)
            save_mesh(expected_case_dir / "D.obj", v_d, f_d)
        if save_expected_result and expected_case_dir is not None:
            save_mesh(expected_case_dir / "Z.obj", v_z, f_z)

        c_path = prep_dir / "C.obj"
        d_path = prep_dir / "D.obj"
        z_path = prep_dir / "Z.obj"
        save_mesh(c_path, v_c, f_c)
        save_mesh(d_path, v_d, f_d)
        save_mesh(z_path, v_z, f_z)

        prepared.update(
            {
                "c_path": str(c_path),
                "d_path": str(d_path),
                "z_path": str(z_path),
                "spheres_a": [
                    {"sphere_id": int(s["sphere_id"]), "center": s["center"].tolist(), "radius": float(s["radius"])}
                    for s in spheres_a
                ],
                "spheres_b": [
                    {"sphere_id": int(s["sphere_id"]), "center": s["center"].tolist(), "radius": float(s["radius"])}
                    for s in spheres_b
                ],
                "cut_stats": {
                    "C": cuts_a,
                    "D": cuts_b,
                    "Z_after_A": z_a_stats,
                    "Z_after_B": z_b_stats,
                },
                "mesh_stats": {
                    "A": mesh_stats(v_a, f_a),
                    "B": mesh_stats(v_b, f_b),
                    "X": mesh_stats(v_x, f_x),
                    "C": mesh_stats(v_c, f_c),
                    "D": mesh_stats(v_d, f_d),
                    "Z": mesh_stats(v_z, f_z),
                },
            }
        )
        if expected_case_dir is not None and debug_requested:
            prepared["expected_results_dir"] = str(expected_case_dir)
    except Exception as exc:
        prepared["status"] = "error"
        prepared["error"] = str(exc)
    return prepared


def prepare_case_data_with_limits(case: Dict[str, Any], config: Dict[str, Any], prep_dir: Path) -> Dict[str, Any]:
    limits_cfg = config.get("limits", {})
    timeout_limit_sec = float(limits_cfg.get("preparation_timeout_sec", 0.0))
    memory_limit_mb = float(limits_cfg.get("preparation_memory_limit_mb", 0.0))
    check_interval = float(limits_cfg.get("memory_check_interval_sec", MEMORY_CHECK_INTERVAL_SEC))
    if timeout_limit_sec <= 0 and memory_limit_mb <= 0:
        prepared = prepare_case_data(case, config, prep_dir)
        prepared["preparation_peak_rss_mb"] = None
        return prepared

    prep_dir.mkdir(parents=True, exist_ok=True)
    case_json = prep_dir / "_prep_case.json"
    cfg_json = prep_dir / "_prep_config.json"
    out_json = prep_dir / "_prep_result.json"
    _write_json(case_json, case)
    _write_json(cfg_json, config)

    worker_script = Path(__file__).resolve().parents[1] / "scripts" / "prepare_case_worker.py"
    argv = [
        sys.executable,
        str(worker_script),
        "--case-json",
        str(case_json),
        "--config-json",
        str(cfg_json),
        "--prep-dir",
        str(prep_dir),
        "--out-json",
        str(out_json),
    ]
    completed, peak_rss_mb, _, terminated_reason = _run_subprocess_with_limits(
        argv=argv,
        cwd=worker_script.parent,
        timeout_limit_sec=timeout_limit_sec,
        memory_limit_mb=memory_limit_mb,
        check_interval_sec=check_interval,
    )
    if terminated_reason is not None:
        return {
            "case_id": str(case.get("case_id", "unknown")),
            "status": "error",
            "error": f"preparation {terminated_reason}",
            "preparation_peak_rss_mb": peak_rss_mb if peak_rss_mb is not None else None,
            "preparation_status_code": int(completed.returncode),
            "preparation_stdout": completed.stdout,
            "preparation_stderr": completed.stderr,
            "preparation_terminated_reason": terminated_reason,
        }
    if completed.returncode != 0:
        return {
            "case_id": str(case.get("case_id", "unknown")),
            "status": "error",
            "error": f"preparation process failed with status code {completed.returncode}",
            "preparation_peak_rss_mb": peak_rss_mb if peak_rss_mb is not None else None,
            "preparation_status_code": int(completed.returncode),
            "preparation_stdout": completed.stdout,
            "preparation_stderr": completed.stderr,
        }
    if not out_json.exists():
        return {
            "case_id": str(case.get("case_id", "unknown")),
            "status": "error",
            "error": "preparation result file not created",
            "preparation_peak_rss_mb": peak_rss_mb if peak_rss_mb is not None else None,
            "preparation_status_code": int(completed.returncode),
            "preparation_stdout": completed.stdout,
            "preparation_stderr": completed.stderr,
        }
    try:
        with out_json.open("r", encoding="utf-8") as handle:
            prepared = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "case_id": str(case.get("case_id", "unknown")),
            "status": "error",
            "error": f"invalid preparation result payload: {exc}",
            "preparation_peak_rss_mb": peak_rss_mb if peak_rss_mb is not None else None,
            "preparation_status_code": int(completed.returncode),
            "preparation_stdout": completed.stdout,
            "preparation_stderr": completed.stderr,
        }
    if isinstance(prepared, dict):
        prepared["preparation_peak_rss_mb"] = peak_rss_mb if peak_rss_mb is not None else None
        prepared["preparation_status_code"] = int(completed.returncode)
        return prepared
    return {
        "case_id": str(case.get("case_id", "unknown")),
        "status": "error",
        "error": "invalid preparation result payload",
        "preparation_peak_rss_mb": peak_rss_mb if peak_rss_mb is not None else None,
        "preparation_status_code": int(completed.returncode),
        "preparation_stdout": completed.stdout,
        "preparation_stderr": completed.stderr,
    }


def run_case_with_prepared(
    case: Dict[str, Any],
    config: Dict[str, Any],
    output_root: Path,
    prepared: Dict[str, Any],
) -> Dict[str, Any]:
    case_id = str(case["case_id"])
    case_dir = output_root / case_id
    case_dir.mkdir(parents=True, exist_ok=True)
    meta: Dict[str, Any] = {
        "case_id": case_id,
        "input_a": case["input_a"],
        "input_b": case["input_b"],
        "operation": _normalize_op(case["operation"]),
        "status": "ok",
    }
    if prepared.get("status") != "ok":
        meta["status"] = "error"
        meta["error"] = f"case preparation failed: {prepared.get('error', 'unknown error')}"
        if prepared.get("preparation_status_code") is not None:
            meta["preparation_status_code"] = prepared.get("preparation_status_code")
        if prepared.get("preparation_peak_rss_mb") is not None:
            meta["preparation_peak_rss_mb"] = prepared.get("preparation_peak_rss_mb")
        _write_json(case_dir / "case_result.json", meta)
        return meta

    try:
        c_path = Path(str(prepared["c_path"]))
        d_path = Path(str(prepared["d_path"]))
        v_z, f_z = (
            (np.asarray(prepared["v_z"], dtype=np.float64), np.asarray(prepared["f_z"], dtype=np.int64))
            if "v_z" in prepared and "f_z" in prepared
            else load_mesh(Path(str(prepared["z_path"])))
        )

        method_cfg = config["method_under_test"]
        method_input_a = c_path
        method_input_b = d_path
        method_output_y = case_dir / "Y.obj"
        v_y, f_y, method_meta = run_method_under_test(
            method_cfg,
            meta["operation"],
            case_dir,
            method_input_a,
            method_input_b,
            method_output_y,
            limits_cfg=config.get("limits", {}),
        )
        meta["method"] = method_meta

        if v_y is not None and f_y is not None:
            save_mesh(method_output_y, v_y, f_y)
            if bool(method_cfg.get("save_output_meshes", False)):
                mesh_dir = output_root / "meshes"
                mesh_dir.mkdir(parents=True, exist_ok=True)
                method_name = str(method_cfg.get("method_name", "method"))
                mesh_name = f"{case_id}_{meta['operation']}_{method_name}.obj"
                save_mesh(mesh_dir / mesh_name, v_y, f_y)
            meta["metrics"] = evaluate_metrics(v_y, f_y, v_z, f_z, config["metrics"])
        else:
            meta["status"] = "error"
            meta["error"] = "method_under_test failed to produce output"

        if not bool(config.get("debug", {}).get("save_cut_meshes", False)) and method_output_y.exists():
            method_output_y.unlink()

        meta["spheres_a"] = prepared.get("spheres_a", [])
        meta["spheres_b"] = prepared.get("spheres_b", [])
        meta["cut_stats"] = prepared.get("cut_stats", {})
        meta["mesh_stats"] = prepared.get("mesh_stats", {})
        meta["preparation_peak_rss_mb"] = prepared.get("preparation_peak_rss_mb")
        if prepared.get("preparation_status_code") is not None:
            meta["preparation_status_code"] = prepared.get("preparation_status_code")
        if "expected_results_dir" in prepared:
            meta["expected_results_dir"] = prepared["expected_results_dir"]
    except Exception as exc:
        meta["status"] = "error"
        meta["error"] = str(exc)

    _write_json(case_dir / "case_result.json", meta)
    return meta


def run_case(case: Dict[str, Any], config: Dict[str, Any], output_root: Path) -> Dict[str, Any]:
    prep_dir = output_root / "_prepared" / str(case["case_id"])
    prepared = prepare_case_data_with_limits(case, config, prep_dir)
    if prepared.get("status") == "ok":
        v_z, f_z = load_mesh(Path(str(prepared["z_path"])))
        prepared["v_z"] = v_z
        prepared["f_z"] = f_z
    return run_case_with_prepared(case, config, output_root, prepared)


def run_manual_case(
    input_a: Path,
    input_b: Path,
    operation: str,
    case_name: str,
    config: Dict[str, Any],
    output_root: Path,
) -> Dict[str, Any]:
    seed = int(config["seed"])
    rng = np.random.default_rng(seed)
    case = {
        "case_id": case_name,
        "input_a": str(input_a),
        "input_b": str(input_b),
        "operation": _normalize_op(operation),
        "seed_spheres_a": int(rng.integers(0, 2**31 - 1)),
        "seed_spheres_b": int(rng.integers(0, 2**31 - 1)),
    }
    _write_json(output_root / "cases_manifest.json", [case])
    return run_case(case, config, output_root)


def run_dataset_benchmark(
    dataset_dir: Path,
    config: Dict[str, Any],
    output_root: Path,
    update_stats: bool = False,
) -> List[Dict[str, Any]]:
    dataset_manifest_dir = output_root / "manifests"
    mesh_records = get_dataset_mesh_stats(dataset_dir, update_stats=update_stats)
    _write_json(dataset_manifest_dir / "mesh_manifest.json", mesh_records)
    _write_csv(dataset_manifest_dir / "mesh_manifest.csv", mesh_records)

    candidate_records = [
        r
        for r in mesh_records
        if r.get("status") == "ok" and bool(r.get("is_closed", False)) and bool(r.get("is_manifold", False))
    ]
    cases = sample_random_cases(
        candidate_records,
        int(config["num_pairs"]),
        list(config["operations"]),
        int(config["seed"]),
    )
    _write_json(output_root / "cases_manifest.json", cases)
    _write_csv(output_root / "cases_manifest.csv", cases)
    input_pairs = [
        {
            "case_id": c["case_id"],
            "input_a": c["input_a"],
            "input_b": c["input_b"],
        }
        for c in cases
    ]
    _write_json(output_root / "input_mesh_pairs.json", input_pairs)
    _write_csv(output_root / "input_mesh_pairs.csv", input_pairs)

    results = [run_case(case, config, output_root) for case in cases]
    return results


def summarize_results(results: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    total = len(results)
    ok = sum(1 for r in results if r.get("status") == "ok")
    failed = total - ok
    summary: Dict[str, Any] = {
        "total_cases": total,
        "ok_cases": ok,
        "failed_cases": failed,
    }
    haus = []
    chamfer = []
    for r in results:
        metrics = r.get("metrics")
        if not metrics:
            continue
        geo = metrics.get("geometry", {})
        h = geo.get("hausdorff")
        c = geo.get("chamfer")
        if isinstance(h, (float, int)) and np.isfinite(h):
            haus.append(float(h))
        if isinstance(c, (float, int)) and np.isfinite(c):
            chamfer.append(float(c))
    if haus:
        summary["hausdorff_mean"] = float(np.mean(haus))
        summary["hausdorff_median"] = float(np.median(haus))
    if chamfer:
        summary["chamfer_mean"] = float(np.mean(chamfer))
        summary["chamfer_median"] = float(np.median(chamfer))
    return summary


def load_config(config_path: Path) -> Dict[str, Any]:
    with config_path.open("r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    cfg.setdefault("seed", 12345)
    cfg.setdefault("num_pairs", 32)
    cfg.setdefault("operations", ["union", "intersection", "difference"])
    cfg.setdefault("spheres", {})
    cfg.setdefault("metrics", {})
    cfg.setdefault("debug", {})
    cfg.setdefault("plots", {})
    cfg.setdefault("limits", {})
    cfg["debug"].setdefault("save_cut_meshes", False)
    cfg["debug"].setdefault("save_expected_result", False)
    cfg["plots"].setdefault("enabled", False)
    cfg["plots"].setdefault("output_subdir", "plots")
    cfg["plots"].setdefault("dpi", 150)
    cfg["plots"].setdefault("formats", ["pdf"])
    cfg["plots"].setdefault("complexity_bins", 80)
    cfg["plots"].setdefault("export_timeout_sec", 30.0)
    cfg["plots"].setdefault("fallback_html_on_export_failure", True)
    cfg["limits"].setdefault("preparation_timeout_sec", 7200)
    cfg["limits"].setdefault("preparation_memory_limit_mb", 32768)
    cfg["limits"].setdefault("method_timeout_sec", 7200)
    cfg["limits"].setdefault("method_memory_limit_mb", 32768)
    cfg["limits"].setdefault("memory_check_interval_sec", 0.01)
    cfg.setdefault("method_under_test", {})
    return cfg


def load_method_configs(methods_dir: Path) -> List[Dict[str, Any]]:
    method_files = sorted(p for p in methods_dir.glob("*.json") if p.is_file())
    if not method_files:
        raise ValueError(f"No method config files found in: {methods_dir}")
    methods: List[Dict[str, Any]] = []
    for method_file in method_files:
        with method_file.open("r", encoding="utf-8") as handle:
            method_cfg = json.load(handle)
        if not isinstance(method_cfg, dict):
            raise ValueError(f"Method config must be a JSON object: {method_file}")
        if "exe_path" not in method_cfg or "operation_map" not in method_cfg or "output_path" not in method_cfg:
            raise ValueError(
                f"Method config missing required fields (exe_path, operation_map, output_path): {method_file}"
            )
        if "save_output_meshes" not in method_cfg:
            method_cfg["save_output_meshes"] = False
        method_cfg["method_name"] = method_file.stem
        method_cfg["config_path"] = str(method_file)
        methods.append(method_cfg)
    return methods


def save_global_outputs(output_root: Path, config: Dict[str, Any], results: Sequence[Dict[str, Any]]) -> None:
    _write_json(output_root / "config_snapshot.json", config)
    _write_json(output_root / "results_cases.json", list(results))
    _write_csv(output_root / "results_cases.csv", list(results))
    summary = summarize_results(results)
    _write_json(output_root / "results_summary.json", summary)
