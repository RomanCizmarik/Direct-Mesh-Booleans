from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None


SUPPORTED_MESH_SUFFIXES = {".obj", ".off", ".stl", ".ply"}


def list_mesh_paths(dataset_dir: Path) -> List[Path]:
    return sorted(
        p for p in dataset_dir.rglob("*") if p.is_file() and p.suffix.lower() in SUPPORTED_MESH_SUFFIXES
    )


def base_record(mesh_path: Path, mesh_id: str) -> Dict[str, Any]:
    size_bytes = int(mesh_path.stat().st_size)
    return {
        "mesh_id": mesh_id,
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
    }


def _to_serializable(obj: Any) -> Any:
    if np is not None:
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(v) for v in obj]
    return obj


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(_to_serializable(payload), handle, indent=2)


def write_csv(path: Path, rows: List[Dict[str, Any]]) -> None:
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


def default_output_paths(
    dataset_dir: Path,
    suffix: str,
    json_path: Optional[Path],
    csv_path: Optional[Path],
) -> Tuple[Path, Path]:
    out_json = json_path.resolve() if json_path is not None else (dataset_dir / f"mesh_stats_{suffix}.json")
    out_csv = csv_path.resolve() if csv_path is not None else (dataset_dir / f"mesh_stats_{suffix}.csv")
    return out_json, out_csv
