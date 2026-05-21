from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
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


SUPPORTED_MESH_SUFFIXES = {".obj", ".off", ".stl", ".ply"}
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


def sample_random_cases(
    mesh_records: Sequence[Dict[str, Any]],
    num_pairs: int,
    operations: Sequence[str],
    seed: int,
) -> List[Dict[str, Any]]:
    valid = [r for r in mesh_records if r.get("status") == "ok"]
    if len(valid) < 2:
        raise ValueError("Need at least 2 valid meshes to sample pairs.")
    if not operations:
        raise ValueError("At least one operation must be provided.")

    rng = np.random.default_rng(seed)
    normalized_ops = [_normalize_op(op) for op in operations]
    cases: List[Dict[str, Any]] = []
    for i in range(num_pairs):
        ia, ib = rng.choice(len(valid), size=2, replace=False)
        op = normalized_ops[int(rng.integers(0, len(normalized_ops)))]
        cases.append(
            {
                "case_id": f"case_{i:05d}",
                "input_a": valid[int(ia)]["path"],
                "input_b": valid[int(ib)]["path"],
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


def run_method_under_test(
    method_cfg: Dict[str, Any],
    operation: str,
    case_dir: Path,
    v_c: np.ndarray,
    f_c: np.ndarray,
    v_d: np.ndarray,
    f_d: np.ndarray,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Dict[str, Any]]:
    op = _normalize_op(operation)
    t0 = time.perf_counter()
    exe_path = method_cfg.get("exe_path")
    operation_map = method_cfg.get("operation_map")
    if not isinstance(exe_path, str) or not exe_path.strip():
        raise ValueError("Method config must provide non-empty 'exe_path'.")
    if not isinstance(operation_map, dict):
        raise ValueError("Method config must provide 'operation_map' dictionary.")

    in_a = case_dir / "C.obj"
    in_b = case_dir / "D.obj"
    out_y = case_dir / "Y.obj"
    if out_y.exists():
        out_y.unlink()
    save_mesh(in_a, v_c, f_c)
    save_mesh(in_b, v_d, f_d)

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
        case_dir / "X.obj",
        case_dir / "Z.obj",
    ]
    attempts: List[Dict[str, Any]] = []

    def run_attempt(argv: List[str], label: str) -> subprocess.CompletedProcess:
        start_attempt = time.time()
        proc = subprocess.run(argv, capture_output=True, text=True, cwd=str(working_directory))
        attempts.append(
            {
                "label": label,
                "argv": argv,
                "status_code": proc.returncode,
                "stdout": proc.stdout,
                "stderr": proc.stderr,
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
    }
    if not attempts:
        meta["status"] = "error"
        meta["error"] = "No process attempts executed."
        return None, None, meta

    last = attempts[-1]
    meta["status_code"] = last["status_code"]
    if last["status_code"] != 0 or not out_y.exists():
        meta["status"] = "error"
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


def run_case(case: Dict[str, Any], config: Dict[str, Any], output_root: Path) -> Dict[str, Any]:
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

    try:
        v_a, f_a = load_mesh(Path(case["input_a"]))
        v_b, f_b = load_mesh(Path(case["input_b"]))
        save_mesh(case_dir / "A.obj", v_a, f_a)
        save_mesh(case_dir / "B.obj", v_b, f_b)

        v_x, f_x, j_x = run_boolean(v_a, f_a, v_b, f_b, meta["operation"])
        save_mesh(case_dir / "X.obj", v_x, f_x)
        np.save(case_dir / "X_birth_indices.npy", j_x)
        source_tag = (j_x >= f_a.shape[0]).astype(np.int64)  # 0=A, 1=B

        sphere_cfg = config["spheres"]
        spheres_a = sample_cut_spheres(v_a, f_a, sphere_cfg, int(case["seed_spheres_a"]))
        spheres_b = sample_cut_spheres(v_b, f_b, sphere_cfg, int(case["seed_spheres_b"]))

        for s in spheres_a:
            save_mesh(case_dir / "cutters" / "A" / f"sphere_{s['sphere_id']:02d}.obj", s["v"], s["f"])
        for s in spheres_b:
            save_mesh(case_dir / "cutters" / "B" / f"sphere_{s['sphere_id']:02d}.obj", s["v"], s["f"])

        prov_a = np.arange(f_a.shape[0], dtype=np.int64)
        v_c, f_c, _, cuts_a = apply_sphere_cuts(
            v_a,
            f_a,
            prov_a,
            spheres_a,
            target_tag=None,
            save_steps_dir=case_dir / "C_cut_steps",
        )
        save_mesh(case_dir / "C.obj", v_c, f_c)

        prov_b = np.arange(f_b.shape[0], dtype=np.int64)
        v_d, f_d, _, cuts_b = apply_sphere_cuts(
            v_b,
            f_b,
            prov_b,
            spheres_b,
            target_tag=None,
            save_steps_dir=case_dir / "D_cut_steps",
        )
        save_mesh(case_dir / "D.obj", v_d, f_d)

        v_z_a, f_z_a, tag_after_a, z_a_stats = apply_sphere_cuts(
            v_x,
            f_x,
            source_tag,
            spheres_a,
            target_tag=0,
            save_steps_dir=case_dir / "Z_cut_A_steps",
        )
        v_z, f_z, tag_final, z_b_stats = apply_sphere_cuts(
            v_z_a,
            f_z_a,
            tag_after_a,
            spheres_b,
            target_tag=1,
            save_steps_dir=case_dir / "Z_cut_B_steps",
        )
        _ = tag_final
        save_mesh(case_dir / "Z.obj", v_z, f_z)

        method_cfg = config["method_under_test"]
        v_y, f_y, method_meta = run_method_under_test(
            method_cfg,
            meta["operation"],
            case_dir,
            v_c,
            f_c,
            v_d,
            f_d,
        )
        meta["method"] = method_meta

        if v_y is not None and f_y is not None:
            save_mesh(case_dir / "Y.obj", v_y, f_y)
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

        meta["spheres_a"] = [
            {"sphere_id": int(s["sphere_id"]), "center": s["center"].tolist(), "radius": float(s["radius"])}
            for s in spheres_a
        ]
        meta["spheres_b"] = [
            {"sphere_id": int(s["sphere_id"]), "center": s["center"].tolist(), "radius": float(s["radius"])}
            for s in spheres_b
        ]
        meta["cut_stats"] = {
            "C": cuts_a,
            "D": cuts_b,
            "Z_after_A": z_a_stats,
            "Z_after_B": z_b_stats,
        }
        meta["mesh_stats"] = {
            "A": mesh_stats(v_a, f_a),
            "B": mesh_stats(v_b, f_b),
            "X": mesh_stats(v_x, f_x),
            "C": mesh_stats(v_c, f_c),
            "D": mesh_stats(v_d, f_d),
            "Z": mesh_stats(v_z, f_z),
        }
    except Exception as exc:
        meta["status"] = "error"
        meta["error"] = str(exc)

    _write_json(case_dir / "case_result.json", meta)
    return meta


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
) -> List[Dict[str, Any]]:
    dataset_manifest_dir = output_root / "manifests"
    mesh_records = index_dataset(dataset_dir)
    _write_json(dataset_manifest_dir / "mesh_manifest.json", mesh_records)
    _write_csv(dataset_manifest_dir / "mesh_manifest.csv", mesh_records)

    cases = sample_random_cases(
        mesh_records,
        int(config["num_pairs"]),
        list(config["operations"]),
        int(config["seed"]),
    )
    _write_json(output_root / "cases_manifest.json", cases)
    _write_csv(output_root / "cases_manifest.csv", cases)

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

