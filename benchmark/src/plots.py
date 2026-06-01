from __future__ import annotations

import csv
import json
import multiprocessing as mp
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import kaleido  # noqa: F401  # required by plotly.write_image


FONT_FAMILY = "Times New Roman"
DISPLAY_NAME_ALIASES = {
    "direct_mesh_booleans": "DMB",
    "direct_mesh_booleans_extension": "DMB_extension",
    "mine": "DMB",
    "mine_extension": "DMB_extension",
}


def _safe_float(value: Any) -> float:
    if isinstance(value, (int, float)):
        out = float(value)
        return out if np.isfinite(out) else float("nan")
    return float("nan")


def _method_display_name(name: str) -> str:
    key = str(name).strip()
    return DISPLAY_NAME_ALIASES.get(key, key)


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
    for method_name_raw, method_dir in _discover_method_dirs(run_root):
        method_name = _method_display_name(method_name_raw)
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
            chamfer_symmetric = _safe_float(geo.get("chamfer_symmetric"))
            combined_size_mb = _safe_float(method_meta.get("input_size_mb"))
            if not np.isfinite(combined_size_mb):
                input_a = case.get("input_a")
                input_b = case.get("input_b")
                if isinstance(input_a, str) and isinstance(input_b, str):
                    try:
                        combined_size_mb = (Path(input_a).stat().st_size + Path(input_b).stat().st_size) / (1024.0 * 1024.0)
                    except OSError:
                        combined_size_mb = float("nan")
            rows.append(
                {
                    "method": method_name,
                    "method_raw": method_name_raw,
                    "case_id": case.get("case_id"),
                    "operation": case.get("operation"),
                    "status": case.get("status"),
                    "success": 1 if case.get("status") == "ok" else 0,
                    "hausdorff": _safe_float(geo.get("hausdorff")),
                    "chamfer_symmetric": chamfer_symmetric,
                    "chamfer_y_to_z": _safe_float(geo.get("chamfer_y_to_z")),
                    "chamfer_z_to_y": _safe_float(geo.get("chamfer_z_to_y")),
                    "d95_y_to_z": _safe_float(geo.get("d95_y_to_z")),
                    "d95_z_to_y": _safe_float(geo.get("d95_z_to_y")),
                    "runtime_sec": _safe_float(method_meta.get("runtime_sec")),
                    "peak_memory_usage": _safe_float(method_meta.get("peak_rss_mb")),
                    "combined_size_mb": combined_size_mb,
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


def _apply_standard_layout(fig: go.Figure) -> go.Figure:
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=24, family=FONT_FAMILY),
        width=1600,
        height=900,
    )
    fig.update_xaxes(gridcolor="lightgrey", zerolinecolor="lightgrey")
    fig.update_yaxes(gridcolor="lightgrey", zerolinecolor="lightgrey")
    return fig


def _image_export_worker(fig_json: Dict[str, Any], out_path: str, scale: float, queue: Any) -> None:
    try:
        local_fig = go.Figure(fig_json)
        local_fig.write_image(out_path, scale=scale)
        queue.put({"ok": True})
    except Exception as exc:  # pragma: no cover - defensive worker fallback
        queue.put({"ok": False, "error": str(exc)})


def _write_image_with_timeout(fig: go.Figure, out_path: Path, scale: float, timeout_sec: float) -> Tuple[bool, Optional[str]]:
    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(target=_image_export_worker, args=(fig.to_plotly_json(), str(out_path), float(scale), queue))
    proc.start()
    proc.join(timeout=float(timeout_sec))
    if proc.is_alive():
        proc.terminate()
        proc.join()
        return False, f"image export timeout after {timeout_sec:.1f}s ({out_path.name})"
    if not queue.empty():
        result = queue.get()
        if bool(result.get("ok", False)):
            return True, None
        return False, str(result.get("error", "image export failed"))
    if proc.exitcode == 0:
        return True, None
    return False, f"image export failed with exit code {proc.exitcode}"


def _save_multi_format(
    fig: go.Figure,
    base_path: Path,
    formats: Sequence[str],
    dpi: int,
    export_timeout_sec: float,
    fallback_html_on_failure: bool,
) -> Tuple[List[str], List[str]]:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    created: List[str] = []
    warnings: List[str] = []
    scale = max(1.0, float(dpi) / 96.0)
    html_fallback_created = False
    for fmt in formats:
        fmt_clean = str(fmt).strip().lower()
        if not fmt_clean:
            continue
        out = base_path.with_suffix(f".{fmt_clean}")
        if fmt_clean in {"html", "htm"}:
            fig.write_html(str(out), include_plotlyjs="cdn")
            created.append(str(out))
            html_fallback_created = True
            continue

        ok, err = _write_image_with_timeout(fig, out, scale=scale, timeout_sec=export_timeout_sec)
        if ok:
            created.append(str(out))
            continue

        warnings.append(f"{base_path.name}.{fmt_clean}: {err}")
        if fallback_html_on_failure and not html_fallback_created:
            html_out = base_path.with_suffix(".html")
            fig.write_html(str(html_out), include_plotlyjs="cdn")
            created.append(str(html_out))
            html_fallback_created = True
    return created, warnings


def _plot_metric_box(rows: Sequence[Dict[str, Any]], metric: str, y_title: str) -> go.Figure:
    fig = go.Figure()
    methods = sorted({str(r["method"]) for r in rows})
    colors = px.colors.qualitative.D3
    for idx, method in enumerate(methods):
        vals = [
            float(r[metric])
            for r in rows
            if r["method"] == method and r["success"] == 1 and np.isfinite(float(r[metric]))
        ]
        if not vals:
            continue
        fig.add_trace(
            go.Box(
                y=vals,
                name=method,
                boxpoints=False,
                marker_color=colors[idx % len(colors)],
            )
        )
    fig.update_layout(
        xaxis_title="Method",
        yaxis_title=y_title,
        legend_title="Method",
    )
    return _apply_standard_layout(fig)


def _quantile_bin_medians(x_vals: np.ndarray, y_vals: np.ndarray, bins: int) -> Tuple[np.ndarray, np.ndarray]:
    if x_vals.size == 0 or y_vals.size == 0:
        return np.zeros((0,), dtype=np.float64), np.zeros((0,), dtype=np.float64)
    bins = max(1, int(bins))
    quantiles = np.linspace(0.0, 1.0, bins + 1)
    edges = np.quantile(x_vals, quantiles)
    edges = np.unique(edges)
    if edges.size <= 1:
        return np.array([float(np.median(x_vals))]), np.array([float(np.median(y_vals))])
    assignments = np.searchsorted(edges, x_vals, side="right") - 1
    assignments = np.clip(assignments, 0, edges.size - 2)
    x_med: List[float] = []
    y_med: List[float] = []
    for b in range(edges.size - 1):
        mask = assignments == b
        if not np.any(mask):
            continue
        x_med.append(float(np.median(x_vals[mask])))
        y_med.append(float(np.median(y_vals[mask])))
    return np.asarray(x_med, dtype=np.float64), np.asarray(y_med, dtype=np.float64)


def _plot_complexity_scatter(
    rows: Sequence[Dict[str, Any]],
    metric: str,
    y_title: str,
    bins: int,
) -> go.Figure:
    fig = go.Figure()
    methods = sorted({str(r["method"]) for r in rows})
    colors = px.colors.qualitative.D3
    for idx, method in enumerate(methods):
        filtered = [
            r
            for r in rows
            if r["method"] == method
            and r["success"] == 1
            and np.isfinite(float(r.get("combined_size_mb", float("nan"))))
            and np.isfinite(float(r.get(metric, float("nan"))))
            and float(r.get("combined_size_mb", 0.0)) > 0.0
        ]
        if not filtered:
            continue
        x = np.asarray([float(r["combined_size_mb"]) for r in filtered], dtype=np.float64)
        y = np.asarray([float(r[metric]) for r in filtered], dtype=np.float64)
        sort_idx = np.argsort(x)
        x = x[sort_idx]
        y = y[sort_idx]

        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode="markers",
                marker=dict(color=colors[idx % len(colors)], size=10, opacity=0.35),
                name=method,
                legendgroup=method,
                showlegend=True,
            )
        )

        x_med, y_med = _quantile_bin_medians(x, y, bins=bins)
        if x_med.size > 0:
            fig.add_trace(
                go.Scatter(
                    x=x_med,
                    y=y_med,
                    mode="lines+markers",
                    marker=dict(color=colors[idx % len(colors)], size=9),
                    line=dict(color=colors[idx % len(colors)], width=3),
                    name=f"{method} median",
                    legendgroup=method,
                    showlegend=False,
                )
            )

    fig.update_layout(
        xaxis_title="Combined input size [MB] (log scale)",
        yaxis_title=y_title,
        legend_title="Method",
    )
    fig.update_xaxes(type="log")
    return _apply_standard_layout(fig)


def generate_standard_plots(run_root: Path, plots_cfg: Dict[str, Any]) -> Dict[str, Any]:
    run_root = run_root.resolve()
    if not run_root.exists():
        raise ValueError(f"Run directory does not exist: {run_root}")

    formats_cfg = plots_cfg.get("formats", ["pdf"])
    formats = formats_cfg if isinstance(formats_cfg, list) else ["pdf"]
    formats = [str(f).strip().lower() for f in formats if str(f).strip()]
    if not formats:
        formats = ["pdf"]
    output_subdir = str(plots_cfg.get("output_subdir", "plots"))
    dpi = int(plots_cfg.get("dpi", 150))
    complexity_bins = int(plots_cfg.get("complexity_bins", 80))
    export_timeout_sec = float(plots_cfg.get("export_timeout_sec", 15.0))
    fallback_html_on_failure = bool(plots_cfg.get("fallback_html_on_export_failure", True))
    plots_dir = run_root / output_subdir
    plots_dir.mkdir(parents=True, exist_ok=True)

    rows = _collect_plot_rows(run_root)
    if not rows:
        raise ValueError(f"No method results found under: {run_root}")
    _write_rows_csv(plots_dir / "plot_data_merged.csv", rows)

    created_files: List[str] = []
    warnings: List[str] = []
    fig = _plot_metric_box(rows, "hausdorff", "Hausdorff distance")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "hausdorff_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_symmetric", "Chamfer distance (symmetric)")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_y_to_z", "Chamfer distance (Y -> Z)")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_y_to_z_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_z_to_y", "Chamfer distance (Z -> Y)")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_z_to_y_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "d95_y_to_z", "Distance D95 (Y -> Z)")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "d95_y_to_z_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "d95_z_to_y", "Distance D95 (Z -> Y)")
    files, warns = _save_multi_format(
        fig,
        plots_dir / "d95_z_to_y_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_complexity_scatter(rows, "runtime_sec", "Elapsed time [s]", bins=complexity_bins)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "runtime_vs_input_size",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_complexity_scatter(rows, "peak_memory_usage", "Peak memory usage [MB]", bins=complexity_bins)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "memory_vs_input_size",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    methods = sorted({str(r["method"]) for r in rows})
    summary = {
        "plots_dir": str(plots_dir),
        "rows": len(rows),
        "methods": methods,
        "created_files": created_files,
        "warnings": warnings,
    }
    with (plots_dir / "plots_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    return summary
