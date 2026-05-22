from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _safe_float(value: Any) -> float:
    if isinstance(value, (int, float)):
        out = float(value)
        return out if np.isfinite(out) else float("nan")
    return float("nan")


def _discover_method_dirs(run_root: Path) -> List[Tuple[str, Path]]:
    summary_path = run_root / "methods_run_summary.json"
    if summary_path.exists():
        with summary_path.open("r", encoding="utf-8") as handle:
            summary = json.load(handle)
        pairs: List[Tuple[str, Path]] = []
        if isinstance(summary, list):
            for item in summary:
                if not isinstance(item, dict):
                    continue
                method = str(item.get("method", "")).strip()
                out_dir_raw = item.get("output_dir")
                if not method or not isinstance(out_dir_raw, str):
                    continue
                out_dir = Path(out_dir_raw)
                if not out_dir.is_absolute():
                    out_dir = run_root / out_dir
                if (out_dir / "results_cases.json").exists():
                    pairs.append((method, out_dir))
        if pairs:
            return pairs

    pairs = []
    for d in sorted(p for p in run_root.iterdir() if p.is_dir()):
        if d.name in {"expected_results", "plots"}:
            continue
        if (d / "results_cases.json").exists():
            pairs.append((d.name, d))
    return pairs


def _collect_plot_rows(run_root: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for method_name, method_dir in _discover_method_dirs(run_root):
        results_path = method_dir / "results_cases.json"
        with results_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if not isinstance(payload, list):
            continue
        for case in payload:
            if not isinstance(case, dict):
                continue
            metrics = case.get("metrics", {}) if isinstance(case.get("metrics"), dict) else {}
            geo = metrics.get("geometry", {}) if isinstance(metrics.get("geometry"), dict) else {}
            method_meta = case.get("method", {}) if isinstance(case.get("method"), dict) else {}
            mesh_stats = case.get("mesh_stats", {}) if isinstance(case.get("mesh_stats"), dict) else {}
            c_stats = mesh_stats.get("C", {}) if isinstance(mesh_stats.get("C"), dict) else {}
            d_stats = mesh_stats.get("D", {}) if isinstance(mesh_stats.get("D"), dict) else {}
            a_stats = mesh_stats.get("A", {}) if isinstance(mesh_stats.get("A"), dict) else {}
            b_stats = mesh_stats.get("B", {}) if isinstance(mesh_stats.get("B"), dict) else {}
            triangles_cd = _safe_float(c_stats.get("num_faces")) + _safe_float(d_stats.get("num_faces"))
            triangles_ab = _safe_float(a_stats.get("num_faces")) + _safe_float(b_stats.get("num_faces"))
            input_triangles = triangles_cd if np.isfinite(triangles_cd) else triangles_ab

            rows.append(
                {
                    "method": method_name,
                    "case_id": case.get("case_id"),
                    "operation": case.get("operation"),
                    "status": case.get("status"),
                    "success": 1 if case.get("status") == "ok" else 0,
                    "hausdorff": _safe_float(geo.get("hausdorff")),
                    "chamfer": _safe_float(geo.get("chamfer")),
                    "runtime_sec": _safe_float(method_meta.get("runtime_sec")),
                    "peak_rss_mb": _safe_float(method_meta.get("peak_rss_mb")),
                    "input_triangles": input_triangles,
                }
            )
    return rows


def _write_rows_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = sorted({k for row in rows for k in row.keys()}) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not headers:
            handle.write("")
            return
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _save_multi_format(fig: plt.Figure, base_path: Path, formats: Sequence[str], dpi: int) -> List[str]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    created: List[str] = []
    for fmt in formats:
        fmt_clean = str(fmt).strip().lower()
        if not fmt_clean:
            continue
        out = base_path.with_suffix(f".{fmt_clean}")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        created.append(str(out))
    plt.close(fig)
    return created


def _plot_metric_boxplot(rows: Sequence[Dict[str, Any]], metric: str, title: str, ylabel: str) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 5))
    methods = sorted({str(r["method"]) for r in rows})
    data: List[List[float]] = []
    labels: List[str] = []
    for method in methods:
        vals = [float(r[metric]) for r in rows if r["method"] == method and np.isfinite(float(r[metric]))]
        if vals:
            data.append(vals)
            labels.append(method)
    if not data:
        ax.text(0.5, 0.5, f"No valid {metric} values found.", ha="center", va="center")
        ax.set_axis_off()
        return fig
    ax.boxplot(data, labels=labels, showfliers=False)
    ax.set_title(title)
    ax.set_xlabel("Method")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3, axis="y")
    return fig


def _plot_complexity_scatter(
    rows: Sequence[Dict[str, Any]],
    metric: str,
    title: str,
    ylabel: str,
) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 5))
    methods = sorted({str(r["method"]) for r in rows})
    colors = plt.cm.get_cmap("tab10", max(len(methods), 1))
    has_points = False
    for idx, method in enumerate(methods):
        x_vals: List[float] = []
        y_vals: List[float] = []
        for row in rows:
            if row["method"] != method:
                continue
            x = _safe_float(row.get("input_triangles"))
            y = _safe_float(row.get(metric))
            if np.isfinite(x) and np.isfinite(y):
                x_vals.append(x)
                y_vals.append(y)
        if x_vals and y_vals:
            has_points = True
            ax.scatter(x_vals, y_vals, alpha=0.75, s=24, color=colors(idx), label=method)
    if not has_points:
        ax.text(0.5, 0.5, f"No valid {metric} values found.", ha="center", va="center")
        ax.set_axis_off()
        return fig
    ax.set_title(title)
    ax.set_xlabel("Input complexity (triangles in C + D)")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.3)
    ax.legend()
    return fig


def generate_standard_plots(run_root: Path, plots_cfg: Dict[str, Any]) -> Dict[str, Any]:
    run_root = run_root.resolve()
    if not run_root.exists():
        raise ValueError(f"Run directory does not exist: {run_root}")

    formats_cfg = plots_cfg.get("formats", ["png"])
    formats = formats_cfg if isinstance(formats_cfg, list) else ["png"]
    formats = [str(f).strip().lower() for f in formats if str(f).strip()]
    if not formats:
        formats = ["png"]
    dpi = int(plots_cfg.get("dpi", 150))
    output_subdir = str(plots_cfg.get("output_subdir", "plots"))
    plots_dir = run_root / output_subdir
    plots_dir.mkdir(parents=True, exist_ok=True)

    rows = _collect_plot_rows(run_root)
    if not rows:
        raise ValueError(f"No method results found under: {run_root}")
    _write_rows_csv(plots_dir / "plot_data_merged.csv", rows)

    created_files: List[str] = []
    fig = _plot_metric_boxplot(rows, "hausdorff", "Hausdorff distance by method", "Hausdorff distance")
    created_files.extend(_save_multi_format(fig, plots_dir / "hausdorff_by_method", formats, dpi))

    fig = _plot_metric_boxplot(rows, "chamfer", "Chamfer distance by method", "Chamfer distance")
    created_files.extend(_save_multi_format(fig, plots_dir / "chamfer_by_method", formats, dpi))

    fig = _plot_complexity_scatter(
        rows,
        "runtime_sec",
        "Time complexity: runtime vs input triangles",
        "Runtime [s]",
    )
    created_files.extend(_save_multi_format(fig, plots_dir / "runtime_complexity", formats, dpi))

    fig = _plot_complexity_scatter(
        rows,
        "peak_rss_mb",
        "Memory complexity: peak RSS vs input triangles",
        "Peak RSS [MB]",
    )
    created_files.extend(_save_multi_format(fig, plots_dir / "memory_complexity", formats, dpi))

    methods = sorted({str(r["method"]) for r in rows})
    summary = {
        "plots_dir": str(plots_dir),
        "rows": len(rows),
        "methods": methods,
        "created_files": created_files,
    }
    with (plots_dir / "plots_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary

