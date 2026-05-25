from __future__ import annotations

import argparse
import gc
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

try:
    import pymeshlab
except ImportError as exc:  # pragma: no cover
    raise RuntimeError(
        "Python package 'pymeshlab' is required for this script. Install it in your environment first."
    ) from exc

from mesh_stats_common import base_record, default_output_paths, list_mesh_paths, write_csv, write_json


def _normalize_key(key: str) -> str:
    return "".join(ch for ch in key.lower() if ch.isalnum())


def _to_bool(value: Any) -> Optional[bool]:
    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        value_l = value.strip().lower()
        if value_l in {"true", "yes", "y", "1"}:
            return True
        if value_l in {"false", "no", "n", "0"}:
            return False
    return None


def _lookup_measure(measures: Dict[str, Any], aliases: List[str]) -> Any:
    if not isinstance(measures, dict):
        return None
    norm_map = {_normalize_key(str(k)): v for k, v in measures.items()}
    for alias in aliases:
        key = _normalize_key(alias)
        if key in norm_map:
            return norm_map[key]
    for alias in aliases:
        key = _normalize_key(alias)
        for candidate, value in norm_map.items():
            if key in candidate:
                return value
    return None


def _run_filter(ms: pymeshlab.MeshSet, name: str) -> Dict[str, Any]:
    if hasattr(ms, name):
        result = getattr(ms, name)()
        if isinstance(result, dict):
            return result
        return {}
    result = ms.apply_filter(name)
    if isinstance(result, dict):
        return result
    return {}


def _compute_record(mesh_path: Path, mesh_id: str) -> Dict[str, Any]:
    record = base_record(mesh_path, mesh_id)
    ms = pymeshlab.MeshSet()
    try:
        ms.load_new_mesh(str(mesh_path))
        mesh = ms.current_mesh()

        num_vertices = int(mesh.vertex_number())
        num_faces = int(mesh.face_number())
        if num_vertices <= 0 or num_faces <= 0:
            record["status"] = "error"
            record["error"] = "empty mesh"
            return record

        bbox = mesh.bounding_box()
        bmin = np.asarray(bbox.min(), dtype=np.float64)
        bmax = np.asarray(bbox.max(), dtype=np.float64)
        dx = float(bmax[0] - bmin[0])
        dy = float(bmax[1] - bmin[1])
        dz = float(bmax[2] - bmin[2])

        topo = _run_filter(ms, "compute_topological_measures")
        geom = _run_filter(ms, "compute_geometric_measures")

        boundary_edges_raw = _lookup_measure(
            topo,
            [
                "boundary_edges",
                "number_boundary_edges",
                "mesh_boundary_edges",
                "boundaryedge",
            ],
        )
        boundary_edges = int(boundary_edges_raw) if boundary_edges_raw is not None else -1

        is_closed_raw = _lookup_measure(topo, ["is_closed", "is_mesh_closed", "watertight"])
        is_closed = _to_bool(is_closed_raw)

        non_manifold_edges_raw = _lookup_measure(
            topo,
            [
                "non_two_manifold_edges",
                "non_manifold_edges",
                "number_non_manifold_edges",
            ],
        )
        non_manifold_vertices_raw = _lookup_measure(
            topo,
            [
                "non_two_manifold_vertices",
                "non_manifold_vertices",
                "number_non_manifold_vertices",
            ],
        )

        is_edge_manifold = bool(int(non_manifold_edges_raw) == 0) if non_manifold_edges_raw is not None else False
        is_vertex_manifold = bool(int(non_manifold_vertices_raw) == 0) if non_manifold_vertices_raw is not None else False

        is_manifold = bool(is_edge_manifold and is_vertex_manifold)
        surface_area_raw = _lookup_measure(geom, ["surface_area", "mesh_surface_area", "area"])
        surface_area = float(surface_area_raw) if surface_area_raw is not None else 0.0

        record.update(
            {
                "num_vertices": num_vertices,
                "num_faces": num_faces,
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "bbox_diag": float(math.sqrt(dx * dx + dy * dy + dz * dz)),
                "boundary_edges": int(boundary_edges),
                "surface_area": surface_area,
                "num_triangles": num_faces,
                "is_closed": bool(is_closed) if is_closed is not None else False,
                "is_manifold": is_manifold,
                "is_edge_manifold": is_edge_manifold,
                "is_vertex_manifold": is_vertex_manifold,
            }
        )
    except Exception as exc:
        record["status"] = "error"
        record["error"] = str(exc)
    finally:
        ms.clear()
    return record


def compute_dataset_mesh_stats(dataset_dir: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for idx, mesh_path in enumerate(list_mesh_paths(dataset_dir)):
        records.append(_compute_record(mesh_path, f"mesh_{idx:05d}"))
        if (idx + 1) % 20 == 0:
            gc.collect()
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute dataset mesh stats using PyMeshLab.")
    parser.add_argument("--dataset-dir", type=Path, required=True, help="Dataset root directory.")
    parser.add_argument(
        "--json-path",
        type=Path,
        default=None,
        help="Output JSON path (default: mesh_stats_pymeshlab.json).",
    )
    parser.add_argument(
        "--csv-path",
        type=Path,
        default=None,
        help="Output CSV path (default: mesh_stats_pymeshlab.csv).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    json_path, csv_path = default_output_paths(dataset_dir, "pymeshlab", args.json_path, args.csv_path)

    records = compute_dataset_mesh_stats(dataset_dir)
    write_json(json_path, records)
    write_csv(csv_path, records)
    print(
        json.dumps(
            {
                "library": "pymeshlab",
                "dataset_dir": str(dataset_dir),
                "mesh_count": len(records),
                "json_path": str(json_path),
                "csv_path": str(csv_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
