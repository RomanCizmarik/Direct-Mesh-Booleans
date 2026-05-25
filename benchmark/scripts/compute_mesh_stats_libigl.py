from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
from typing import Any, Dict, List

import igl
import numpy as np

from mesh_stats_common import base_record, default_output_paths, list_mesh_paths, write_csv, write_json


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


def _compute_record(mesh_path: Path, mesh_id: str) -> Dict[str, Any]:
    print(mesh_path)
    record = base_record(mesh_path, mesh_id)
    try:
        v, f = igl.read_triangle_mesh(str(mesh_path))
        v, I, J, f = igl.remove_duplicate_vertices(v, f, 0)

        if v.size == 0 and f.size == 0:
            record["status"] = "error"
            record["error"] = "empty mesh"
            return record
        if v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or f.shape[1] != 3:
            record["status"] = "error"
            record["error"] = "non-triangular mesh"
            return record

        boundary_raw = igl.boundary_facets(f)
        boundary = boundary_raw[0] if isinstance(boundary_raw, tuple) else boundary_raw
        boundary = np.asarray(boundary, dtype=np.int64)

        bmin = v.min(axis=0)
        bmax = v.max(axis=0)
        area = float(np.asarray(igl.doublearea(v, f), dtype=np.float64).sum() * 0.5)
        is_edge_manifold = _edge_manifold_bool(f)
        is_vertex_manifold = _vertex_manifold_bool(f)
        is_manifold = bool(is_edge_manifold and is_vertex_manifold)
        is_closed = bool(boundary.shape[0] == 0)

        record.update(
            {
                "num_vertices": int(v.shape[0]),
                "num_faces": int(f.shape[0]),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "bbox_diag": float(np.linalg.norm(bmax - bmin)),
                "boundary_edges": int(boundary.shape[0]),
                "surface_area": area,
                "num_triangles": int(f.shape[0]),
                "is_closed": is_closed,
                "is_manifold": is_manifold,
                "is_edge_manifold": bool(is_edge_manifold),
                "is_vertex_manifold": bool(is_vertex_manifold),
            }
        )
    except Exception as exc:
        record["status"] = "error"
        record["error"] = str(exc)
    return record


def compute_dataset_mesh_stats(dataset_dir: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for idx, mesh_path in enumerate(list_mesh_paths(dataset_dir)):
        records.append(_compute_record(mesh_path, f"mesh_{idx:05d}"))
        if (idx + 1) % 20 == 0:
            gc.collect()
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute dataset mesh stats using libigl.")
    parser.add_argument("--dataset-dir", type=Path, required=True, help="Dataset root directory.")
    parser.add_argument("--json-path", type=Path, default=None, help="Output JSON path (default: mesh_stats_libigl.json).")
    parser.add_argument("--csv-path", type=Path, default=None, help="Output CSV path (default: mesh_stats_libigl.csv).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    json_path, csv_path = default_output_paths(dataset_dir, "libigl", args.json_path, args.csv_path)

    records = compute_dataset_mesh_stats(dataset_dir)
    write_json(json_path, records)
    write_csv(csv_path, records)
    print(
        json.dumps(
            {
                "library": "libigl",
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
