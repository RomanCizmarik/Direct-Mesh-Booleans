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
MIN_POS_FLOAT = float(np.nextafter(0.0, 1.0))
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
                    "hausdorff_y_to_z": _safe_float(geo.get("hausdorff_y_to_z")),
                    "hausdorff_z_to_y": _safe_float(geo.get("hausdorff_z_to_y")),
                    "chamfer_symmetric": chamfer_symmetric,
                    "chamfer_y_to_z": _safe_float(geo.get("chamfer_y_to_z")),
                    "chamfer_z_to_y": _safe_float(geo.get("chamfer_z_to_y")),
                    "chamfer_winsorized_symmetric": _safe_float(geo.get("chamfer_winsorized_symmetric")),
                    "chamfer_winsorized_y_to_z": _safe_float(geo.get("chamfer_winsorized_y_to_z")),
                    "chamfer_winsorized_z_to_y": _safe_float(geo.get("chamfer_winsorized_z_to_y")),
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
    static_image_fig: Optional[go.Figure] = None,
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

        fig_for_static = static_image_fig if static_image_fig is not None else fig
        ok, err = _write_image_with_timeout(fig_for_static, out, scale=scale, timeout_sec=export_timeout_sec)
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


def _format_stat_value(value: float) -> str:
    v = float(value)
    av = abs(v)
    if av == 0.0:
        return "0"
    if av < 1e-3 or av >= 1e3:
        return f"{v:.2e}"
    return f"{v:.6g}"


def _spread_three_values(min_v: float, med_v: float, max_v: float, gap: float) -> Tuple[float, float, float]:
    y0 = float(min_v)
    y1 = float(med_v)
    y2 = float(max_v)

    if y1 < y0 + gap:
        y1 = y0 + gap
    if y2 < y1 + gap:
        y2 = y1 + gap

    raw_center = (float(min_v) + float(med_v) + float(max_v)) / 3.0
    adj_center = (y0 + y1 + y2) / 3.0
    shift = raw_center - adj_center
    y0 += shift
    y1 += shift
    y2 += shift

    if y1 < y0 + gap:
        y1 = y0 + gap
    if y2 < y1 + gap:
        y2 = y1 + gap
    return y0, y1, y2


def _spread_three_values_log(min_v: float, med_v: float, max_v: float, gap_log: float) -> Tuple[float, float, float]:
    tiny = MIN_POS_FLOAT
    l0 = float(np.log10(max(min_v, tiny)))
    l1 = float(np.log10(max(med_v, tiny)))
    l2 = float(np.log10(max(max_v, tiny)))

    if l1 < l0 + gap_log:
        l1 = l0 + gap_log
    if l2 < l1 + gap_log:
        l2 = l1 + gap_log

    raw_center = (float(np.log10(max(min_v, tiny))) + float(np.log10(max(med_v, tiny))) + float(np.log10(max(max_v, tiny)))) / 3.0
    adj_center = (l0 + l1 + l2) / 3.0
    shift = raw_center - adj_center
    l0 += shift
    l1 += shift
    l2 += shift

    if l1 < l0 + gap_log:
        l1 = l0 + gap_log
    if l2 < l1 + gap_log:
        l2 = l1 + gap_log
    return float(10.0**l0), float(10.0**l1), float(10.0**l2)


def _add_box_summary_annotations(fig: go.Figure, stats_rows: Sequence[Dict[str, Any]], log_scale: bool = False) -> None:
    if not stats_rows:
        return

    global_min = min(float(s["min"]) for s in stats_rows)
    global_max = max(float(s["max"]) for s in stats_rows)
    if log_scale:
        tiny = MIN_POS_FLOAT
        global_min = max(global_min, tiny)
        global_max = max(global_max, global_min * 1.001)
        span_log = max(np.log10(global_max) - np.log10(global_min), 0.25)
        min_gap_log = max(span_log * 0.07, 0.08)
    else:
        span = max(global_max - global_min, max(abs(global_min), abs(global_max)) * 0.05, 1e-12)
        min_gap = max(span * 0.035, 1e-12)

    adj_min = float("inf")
    adj_max = -float("inf")
    for s in stats_rows:
        method = str(s["method"])
        min_raw = float(s["min"])
        med_raw = float(s["median"])
        max_raw = float(s["max"])
        color = str(s.get("color", "rgba(34, 102, 170, 0.95)"))
        if log_scale:
            min_y, med_y, max_y = _spread_three_values_log(min_raw, med_raw, max_raw, min_gap_log)
        else:
            min_y, med_y, max_y = _spread_three_values(min_raw, med_raw, max_raw, min_gap)
        adj_min = min(adj_min, min_y)
        adj_max = max(adj_max, max_y)

        for label, raw_v, y in (
            ("min", min_raw, min_y),
            ("median", med_raw, med_y),
            ("max", max_raw, max_y),
        ):
            fig.add_annotation(
                x=method,
                y=float(y),
                xref="x",
                yref="y",
                text=f"{label}: {_format_stat_value(raw_v)}",
                showarrow=False,
                xanchor="left",
                yanchor="middle",
                xshift=96,
                bgcolor=color,
                bordercolor=color,
                font=dict(size=14, color="white", family=FONT_FAMILY),
                align="left",
            )

    if log_scale:
        tiny = MIN_POS_FLOAT
        lo = min(max(global_min, tiny), max(adj_min, tiny))
        hi = max(global_max, adj_max)
        lo_log = np.log10(lo) - (span_log * 0.10 + min_gap_log)
        hi_log = np.log10(hi) + (span_log * 0.12 + min_gap_log)
        fig.update_yaxes(range=[float(lo_log), float(hi_log)])
    else:
        y_pad = span * 0.08 + min_gap
        fig.update_yaxes(range=[min(global_min, adj_min) - y_pad, max(global_max, adj_max) + y_pad])
    fig.update_layout(margin=dict(r=280))


def _plot_metric_box(
    rows: Sequence[Dict[str, Any]],
    metric: str,
    y_title: str,
    show_summary_stats: bool = False,
    log_scale: bool = False,
    log_epsilon: float = MIN_POS_FLOAT,
) -> go.Figure:
    fig = go.Figure()
    methods = sorted({str(r["method"]) for r in rows})
    colors = px.colors.qualitative.D3
    stats_rows: List[Dict[str, Any]] = []
    for idx, method in enumerate(methods):
        vals: List[float] = []
        for r in rows:
            if r["method"] != method or r["success"] != 1:
                continue
            val = float(r.get(metric, float("nan")))
            if not np.isfinite(val):
                continue
            if log_scale:
                vals.append(float(max(val, float(log_epsilon))))
            else:
                vals.append(val)
        if not vals:
            continue
        arr = np.asarray(vals, dtype=np.float64)
        stats_rows.append(
            {
                "method": method,
                "min": float(np.min(arr)),
                "median": float(np.median(arr)),
                "max": float(np.max(arr)),
                "color": colors[idx % len(colors)],
            }
        )
        fig.add_trace(
            go.Box(
                y=vals,
                name=method,
                boxpoints=False,
                marker_color=colors[idx % len(colors)],
            )
        )

    if log_scale:
        fig.update_yaxes(type="log")
    if show_summary_stats:
        _add_box_summary_annotations(fig, stats_rows, log_scale=log_scale)

    fig.update_layout(
        xaxis_title="Method",
        yaxis_title=y_title,
        legend_title="Method",
    )
    return _apply_standard_layout(fig)


def _plot_asymmetry_box(
    rows: Sequence[Dict[str, Any]],
    y_to_z_key: str,
    z_to_y_key: str,
    y_title: str,
    eps: float = 1e-12,
    show_summary_stats: bool = False,
) -> go.Figure:
    fig = go.Figure()
    methods = sorted({str(r["method"]) for r in rows})
    colors = px.colors.qualitative.D3
    stats_rows: List[Dict[str, Any]] = []
    for idx, method in enumerate(methods):
        vals: List[float] = []
        for r in rows:
            if r["method"] != method or r["success"] != 1:
                continue
            yz = float(r.get(y_to_z_key, float("nan")))
            zy = float(r.get(z_to_y_key, float("nan")))
            if not np.isfinite(yz) or not np.isfinite(zy):
                continue
            vals.append(float(np.log10((yz + eps) / (zy + eps))))
        if not vals:
            continue
        arr = np.asarray(vals, dtype=np.float64)
        stats_rows.append(
            {
                "method": method,
                "min": float(np.min(arr)),
                "median": float(np.median(arr)),
                "max": float(np.max(arr)),
                "color": colors[idx % len(colors)],
            }
        )
        fig.add_trace(
            go.Box(
                y=vals,
                name=method,
                boxpoints=False,
                marker_color=colors[idx % len(colors)],
            )
        )
    if show_summary_stats:
        _add_box_summary_annotations(fig, stats_rows)
    fig.update_layout(
        xaxis_title="Method",
        yaxis_title=y_title,
        legend_title="Method",
    )
    return _apply_standard_layout(fig)


def _plot_success_rate(rows: Sequence[Dict[str, Any]]) -> go.Figure:
    methods = sorted({str(r["method"]) for r in rows})
    rates: List[float] = []
    labels: List[str] = []
    for method in methods:
        subset = [r for r in rows if r["method"] == method]
        total = len(subset)
        succ = sum(int(r.get("success", 0)) for r in subset)
        rate = (100.0 * succ / total) if total > 0 else 0.0
        rates.append(rate)
        labels.append(f"{succ}/{total}")

    fig = go.Figure(
        data=[
            go.Bar(
                x=methods,
                y=rates,
                text=labels,
                textposition="outside",
                marker_color=px.colors.qualitative.D3[: len(methods)],
            )
        ]
    )
    fig.update_layout(
        xaxis_title="Method",
        yaxis_title="Success rate [%]",
        showlegend=False,
    )
    fig.update_yaxes(range=[0.0, 105.0])
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
    fig = _plot_metric_box(rows, "hausdorff", "Hausdorff distance", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "hausdorff", "Hausdorff distance", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "hausdorff_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "hausdorff_y_to_z", "Hausdorff distance (Y -> Z)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "hausdorff_y_to_z", "Hausdorff distance (Y -> Z)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "hausdorff_y_to_z_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "hausdorff_z_to_y", "Hausdorff distance (Z -> Y)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "hausdorff_z_to_y", "Hausdorff distance (Z -> Y)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "hausdorff_z_to_y_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_symmetric", "Chamfer distance (symmetric)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "chamfer_symmetric", "Chamfer distance (symmetric)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(
        rows,
        "chamfer_symmetric",
        "Chamfer distance (symmetric, log scale)",
        show_summary_stats=False,
        log_scale=True,
    )
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_by_method_log",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(
        rows,
        "chamfer_winsorized_symmetric",
        "Chamfer distance (winsorized symmetric)",
        show_summary_stats=False,
    )
    fig_static = _plot_metric_box(
        rows,
        "chamfer_winsorized_symmetric",
        "Chamfer distance (winsorized symmetric)",
        show_summary_stats=True,
    )
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_winsorized_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_y_to_z", "Chamfer distance (Y -> Z)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "chamfer_y_to_z", "Chamfer distance (Y -> Z)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_y_to_z_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "chamfer_z_to_y", "Chamfer distance (Z -> Y)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "chamfer_z_to_y", "Chamfer distance (Z -> Y)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_z_to_y_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "d95_y_to_z", "Distance D95 (Y -> Z)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "d95_y_to_z", "Distance D95 (Y -> Z)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "d95_y_to_z_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_metric_box(rows, "d95_z_to_y", "Distance D95 (Z -> Y)", show_summary_stats=False)
    fig_static = _plot_metric_box(rows, "d95_z_to_y", "Distance D95 (Z -> Y)", show_summary_stats=True)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "d95_z_to_y_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_asymmetry_box(
        rows,
        "chamfer_y_to_z",
        "chamfer_z_to_y",
        "Chamfer asymmetry log10((Y -> Z + eps)/(Z -> Y + eps))",
        eps=1e-12,
        show_summary_stats=False,
    )
    fig_static = _plot_asymmetry_box(
        rows,
        "chamfer_y_to_z",
        "chamfer_z_to_y",
        "Chamfer asymmetry log10((Y -> Z + eps)/(Z -> Y + eps))",
        eps=1e-12,
        show_summary_stats=True,
    )
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_asymmetry_logratio_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_asymmetry_box(
        rows,
        "chamfer_winsorized_y_to_z",
        "chamfer_winsorized_z_to_y",
        "Winsorized Chamfer asymmetry log10((Y -> Z + eps)/(Z -> Y + eps))",
        eps=1e-12,
        show_summary_stats=False,
    )
    fig_static = _plot_asymmetry_box(
        rows,
        "chamfer_winsorized_y_to_z",
        "chamfer_winsorized_z_to_y",
        "Winsorized Chamfer asymmetry log10((Y -> Z + eps)/(Z -> Y + eps))",
        eps=1e-12,
        show_summary_stats=True,
    )
    files, warns = _save_multi_format(
        fig,
        plots_dir / "chamfer_winsorized_asymmetry_logratio_by_method",
        formats,
        dpi,
        export_timeout_sec,
        fallback_html_on_failure,
        static_image_fig=fig_static,
    )
    created_files.extend(files)
    warnings.extend(warns)

    fig = _plot_success_rate(rows)
    files, warns = _save_multi_format(
        fig,
        plots_dir / "success_rate_by_method",
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
