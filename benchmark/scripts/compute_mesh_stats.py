from __future__ import annotations

import argparse
import gc
import json
import math
from pathlib import Path
from typing import Any, Dict, List

import trimesh

from mesh_stats_common import base_record, list_mesh_paths, write_csv, write_json


def _load_trimesh(mesh_path: Path) -> trimesh.Trimesh:
    loaded = trimesh.load_mesh(str(mesh_path), process=False)
    if isinstance(loaded, trimesh.Scene):
        if not loaded.geometry:
            raise ValueError("empty scene")
        geom = [g for g in loaded.geometry.values() if isinstance(g, trimesh.Trimesh)]
        if not geom:
            raise ValueError("scene contains no mesh geometry")
        loaded = trimesh.util.concatenate(geom)
    if not isinstance(loaded, trimesh.Trimesh):
        raise ValueError(f"Unsupported trimesh payload type: {type(loaded).__name__}")
    if loaded.vertices is None or loaded.faces is None:
        raise ValueError("invalid mesh payload")
    if loaded.faces.ndim != 2 or loaded.faces.shape[1] != 3:
        raise ValueError("non-triangular mesh")
    return loaded


def _compute_record(mesh_path: Path, mesh_id: str) -> Dict[str, Any]:
    record = base_record(mesh_path, mesh_id)
    mesh: trimesh.Trimesh | None = None
    try:
        mesh = _load_trimesh(mesh_path)
        mesh.merge_vertices()

        edge_groups = trimesh.grouping.group_rows(mesh.edges_sorted, require_count=None)
        boundary_edges = int(sum(1 for g in edge_groups if len(g) == 1))
        is_closed = bool(boundary_edges == 0)
        is_edge_manifold = bool(all(len(g) <= 2 for g in edge_groups))
        is_vertex_manifold = bool(mesh.is_winding_consistent)
        is_manifold = bool(is_edge_manifold and is_vertex_manifold)

        bounds = mesh.bounds
        bmin = bounds[0]
        bmax = bounds[1]
        dx = float(bmax[0] - bmin[0])
        dy = float(bmax[1] - bmin[1])
        dz = float(bmax[2] - bmin[2])

        record.update(
            {
                "num_vertices": int(mesh.vertices.shape[0]),
                "num_faces": int(mesh.faces.shape[0]),
                "bbox_min": bmin.tolist(),
                "bbox_max": bmax.tolist(),
                "bbox_diag": float(math.sqrt(dx * dx + dy * dy + dz * dz)),
                "boundary_edges": boundary_edges,
                "surface_area": float(mesh.area),
                "num_triangles": int(mesh.faces.shape[0]),
                "is_closed": is_closed,
                "is_manifold": is_manifold,
                "is_edge_manifold": is_edge_manifold,
                "is_vertex_manifold": is_vertex_manifold,
            }
        )
    except Exception as exc:
        record["status"] = "error"
        record["error"] = str(exc)
    finally:
        if mesh is not None:
            mesh._cache.clear()
    return record


def compute_dataset_mesh_stats(dataset_dir: Path) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for idx, mesh_path in enumerate(list_mesh_paths(dataset_dir)):
        records.append(_compute_record(mesh_path, f"mesh_{idx:05d}"))
        if (idx + 1) % 20 == 0:
            gc.collect()
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute dataset mesh stats using trimesh.")
    parser.add_argument("--dataset-dir", type=Path, required=True, help="Dataset root directory.")
    parser.add_argument("--json-path", type=Path, default=None, help="Output JSON path (default: mesh_stats.json).")
    parser.add_argument("--csv-path", type=Path, default=None, help="Output CSV path (default: mesh_stats.csv).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    dataset_dir = args.dataset_dir.resolve()
    json_path = args.json_path.resolve() if args.json_path is not None else (dataset_dir / "mesh_stats.json")
    csv_path = args.csv_path.resolve() if args.csv_path is not None else (dataset_dir / "mesh_stats.csv")

    records = compute_dataset_mesh_stats(dataset_dir)
    write_json(json_path, records)
    write_csv(csv_path, records)
    print(
        json.dumps(
            {
                "library": "trimesh",
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
