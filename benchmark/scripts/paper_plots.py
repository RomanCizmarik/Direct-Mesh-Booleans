from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


REPO_BENCHMARK_DIR = Path(__file__).resolve().parents[1]
DEFAULT_RUN_A = REPO_BENCHMARK_DIR / "artifacts" / "10k_1000_pairs_2_3_10_15_v3"
DEFAULT_RUN_B = REPO_BENCHMARK_DIR / "artifacts" / "10k_1000_pairs_10_20_3_5_v3"
DEFAULT_OUT = REPO_BENCHMARK_DIR / "artifacts" / "paper_plots"
DEFAULT_LOG_FLOOR = 1e-6

FONT_FAMILY = "Times New Roman"
MIN_POS_FLOAT = float(np.nextafter(0.0, 1.0))
DISPLAY_NAME_ALIASES = {
    "direct_mesh_booleans": "DMB",
    "direct_mesh_booleans_extension": "DMB_extension",
    "mine": "DMB",
    "mine_extension": "DMB_extension",
}


def _safe_float(value: Any) -> float:
    try:
        out = float(value)
        return out if np.isfinite(out) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _safe_int(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _method_display_name(name: str) -> str:
    key = str(name).strip()
    return DISPLAY_NAME_ALIASES.get(key, key)


def _regime_label_from_config(run_dir: Path) -> str:
    cfg_path = run_dir / "main_config_effective.json"
    if not cfg_path.exists():
        return run_dir.name
    try:
        cfg = json.loads(cfg_path.read_text(encoding="utf-8"))
        spheres = cfg.get("spheres", {})
        cmin = int(spheres.get("count_min"))
        cmax = int(spheres.get("count_max"))
        rmin = float(spheres.get("radius_min_frac"))
        rmax = float(spheres.get("radius_max_frac"))
        return f"{cmin}-{cmax} spheres, r={100*rmin:.0f}-{100*rmax:.0f}%"
    except Exception:
        return run_dir.name


def _load_plot_rows(run_dir: Path) -> List[Dict[str, Any]]:
    csv_path = run_dir / "plots" / "plot_data_merged.csv"
    if not csv_path.exists():
        raise ValueError(f"Missing plot_data_merged.csv: {csv_path}")
    rows: List[Dict[str, Any]] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            method_raw = str(row.get("method_raw", row.get("method", ""))).strip()
            method = _method_display_name(str(row.get("method", method_raw)))
            success = _safe_int(row.get("success"))
            if success not in (0, 1):
                success = 1 if str(row.get("status", "")).lower() == "ok" else 0
            rows.append(
                {
                    "method": method,
                    "method_raw": method_raw,
                    "success": success,
                    "chamfer_symmetric": _safe_float(row.get("chamfer_symmetric")),
                    "chamfer_y_to_z": _safe_float(row.get("chamfer_y_to_z")),
                    "chamfer_z_to_y": _safe_float(row.get("chamfer_z_to_y")),
                    "hausdorff": _safe_float(row.get("hausdorff")),
                }
            )
    return rows


def _filter_metric_values(
    rows: Sequence[Dict[str, Any]],
    method: str,
    metric: str,
    log_axis: bool = False,
    log_floor: float = MIN_POS_FLOAT,
) -> np.ndarray:
    vals = [
        float(r.get(metric, float("nan")))
        for r in rows
        if str(r.get("method", "")) == method and int(r.get("success", 0)) == 1 and np.isfinite(float(r.get(metric, float("nan"))))
    ]
    if not vals:
        return np.zeros((0,), dtype=np.float64)
    arr = np.asarray(vals, dtype=np.float64)
    if log_axis:
        arr = np.maximum(arr, float(log_floor))
    return arr


def _apply_layout(fig: go.Figure) -> go.Figure:
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=20, family=FONT_FAMILY),
        width=1800,
        height=850,
        margin=dict(l=80, r=40, t=90, b=80),
    )
    fig.update_xaxes(gridcolor="lightgrey", zerolinecolor="lightgrey")
    fig.update_yaxes(gridcolor="lightgrey", zerolinecolor="lightgrey")
    return fig


def _hex_to_rgb(color: str) -> Tuple[int, int, int]:
    value = str(color).strip()
    if value.startswith("#") and len(value) == 7:
        return int(value[1:3], 16), int(value[3:5], 16), int(value[5:7], 16)
    return 127, 127, 127


def _blend_with_white(color: str, ratio: float) -> str:
    ratio = min(max(float(ratio), 0.0), 1.0)
    r, g, b = _hex_to_rgb(color)
    rr = int(round(r + (255 - r) * ratio))
    gg = int(round(g + (255 - g) * ratio))
    bb = int(round(b + (255 - b) * ratio))
    return f"rgb({rr},{gg},{bb})"


def _blend_with_black(color: str, ratio: float) -> str:
    ratio = min(max(float(ratio), 0.0), 1.0)
    r, g, b = _hex_to_rgb(color)
    rr = int(round(r * (1.0 - ratio)))
    gg = int(round(g * (1.0 - ratio)))
    bb = int(round(b * (1.0 - ratio)))
    return f"rgb({rr},{gg},{bb})"


def _save_figure(fig: go.Figure, out_base: Path, formats: Sequence[str], dpi: int) -> List[str]:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    created: List[str] = []
    scale = max(1.0, float(dpi) / 96.0)
    for fmt in formats:
        fmt_clean = str(fmt).strip().lower()
        if not fmt_clean:
            continue
        out_path = out_base.with_suffix(f".{fmt_clean}")
        if fmt_clean in {"html", "htm"}:
            fig.write_html(str(out_path), include_plotlyjs="cdn")
        else:
            fig.write_image(str(out_path), scale=scale)
        created.append(str(out_path))
    return created


def _build_ecdf_figure(
    run_a_rows: Sequence[Dict[str, Any]],
    run_b_rows: Sequence[Dict[str, Any]],
    run_a_label: str,
    run_b_label: str,
    methods: Sequence[str],
    metric: str,
    metric_title: str,
    log_x: bool,
    colors: Dict[str, str],
    log_floor: float = MIN_POS_FLOAT,
) -> go.Figure:
    fig = make_subplots(rows=1, cols=2, subplot_titles=[run_a_label, run_b_label], shared_yaxes=True)
    for col_idx, rows in enumerate((run_a_rows, run_b_rows), start=1):
        for method in methods:
            vals = _filter_metric_values(rows, method, metric, log_axis=log_x, log_floor=log_floor)
            if vals.size == 0:
                continue
            vals = np.sort(vals)
            y = np.arange(1, vals.size + 1, dtype=np.float64) / float(vals.size)
            fig.add_trace(
                go.Scatter(
                    x=vals,
                    y=y,
                    mode="lines",
                    name=method,
                    legendgroup=method,
                    showlegend=(col_idx == 1),
                    line=dict(color=colors[method], width=3),
                ),
                row=1,
                col=col_idx,
            )
        fig.update_xaxes(title_text=f"{metric_title}", row=1, col=col_idx)
        if log_x:
            fig.update_xaxes(type="log", row=1, col=col_idx)
    fig.update_yaxes(title_text="Fraction of successful cases", row=1, col=1, range=[0.0, 1.0])
    fig.update_layout(legend_title="Method")
    return _apply_layout(fig)


def _build_grouped_chamfer_ecdf_figure(
    run_a_rows: Sequence[Dict[str, Any]],
    run_b_rows: Sequence[Dict[str, Any]],
    run_a_label: str,
    run_b_label: str,
    methods: Sequence[str],
    colors: Dict[str, str],
    log_x: bool = False,
    log_floor: float = MIN_POS_FLOAT,
) -> go.Figure:
    metric_specs: Sequence[Tuple[str, str]] = [
        ("chamfer_symmetric", "Chamfer distance (symmetric)"),
        ("chamfer_y_to_z", "Chamfer distance (Y \u2192 Z)"),
        ("chamfer_z_to_y", "Chamfer distance (Z \u2192 Y)"),
    ]
    subplot_titles: List[str] = []
    for _, metric_title in metric_specs:
        subplot_titles.extend([f"{metric_title} \u2014 {run_a_label}", f"{metric_title} \u2014 {run_b_label}"])

    fig = make_subplots(
        rows=len(metric_specs),
        cols=2,
        subplot_titles=subplot_titles,
        shared_yaxes=True,
        vertical_spacing=0.08,
        horizontal_spacing=0.06,
    )
    run_cols = [(1, run_a_rows), (2, run_b_rows)]
    for row_idx, (metric, metric_title) in enumerate(metric_specs, start=1):
        for col_idx, rows in run_cols:
            for method in methods:
                vals = _filter_metric_values(rows, method, metric, log_axis=log_x, log_floor=log_floor)
                if vals.size == 0:
                    continue
                vals = np.sort(vals)
                y = np.arange(1, vals.size + 1, dtype=np.float64) / float(vals.size)
                fig.add_trace(
                    go.Scatter(
                        x=vals,
                        y=y,
                        mode="lines",
                        name=method,
                        legendgroup=method,
                        showlegend=(row_idx == 1 and col_idx == 1),
                        line=dict(color=colors[method], width=2.5),
                    ),
                    row=row_idx,
                    col=col_idx,
                )
            fig.update_xaxes(title_text=metric_title, row=row_idx, col=col_idx)
            if log_x:
                fig.update_xaxes(type="log", row=row_idx, col=col_idx)
        fig.update_yaxes(title_text="Fraction of successful cases", row=row_idx, col=1, range=[0.0, 1.0])
    fig.update_layout(legend_title="Method")
    fig = _apply_layout(fig)
    fig.update_layout(height=max(1300, 420 * len(metric_specs)), width=1900, margin=dict(l=80, r=40, t=110, b=70))
    return fig


def _build_log_box_figure(
    run_a_rows: Sequence[Dict[str, Any]],
    run_b_rows: Sequence[Dict[str, Any]],
    run_a_label: str,
    run_b_label: str,
    methods: Sequence[str],
    metric: str,
    metric_title: str,
    colors: Dict[str, str],
    log_floor: float = MIN_POS_FLOAT,
) -> go.Figure:
    fig = make_subplots(rows=1, cols=2, subplot_titles=[run_a_label, run_b_label], shared_yaxes=True)
    for col_idx, rows in enumerate((run_a_rows, run_b_rows), start=1):
        for method in methods:
            vals = _filter_metric_values(rows, method, metric, log_axis=True, log_floor=log_floor)
            if vals.size == 0:
                continue
            fig.add_trace(
                go.Box(
                    y=vals,
                    name=method,
                    legendgroup=method,
                    showlegend=(col_idx == 1),
                    boxpoints=False,
                    marker_color=colors[method],
                ),
                row=1,
                col=col_idx,
            )
        fig.update_xaxes(title_text="Method", row=1, col=col_idx)
        fig.update_yaxes(type="log", row=1, col=col_idx)
    fig.update_yaxes(title_text=f"{metric_title} (log scale)", row=1, col=1)
    fig.update_layout(legend_title="Method")
    return _apply_layout(fig)


def _build_violin_figure(
    run_a_rows: Sequence[Dict[str, Any]],
    run_b_rows: Sequence[Dict[str, Any]],
    run_a_label: str,
    run_b_label: str,
    methods: Sequence[str],
    metric: str,
    metric_title: str,
    colors: Dict[str, str],
    log_y: bool = False,
    log_floor: float = MIN_POS_FLOAT,
) -> go.Figure:
    fig = make_subplots(rows=1, cols=2, subplot_titles=[run_a_label, run_b_label], shared_yaxes=True)
    for col_idx, rows in enumerate((run_a_rows, run_b_rows), start=1):
        for method in methods:
            vals = _filter_metric_values(rows, method, metric, log_axis=log_y, log_floor=log_floor)
            if vals.size == 0:
                continue
            fig.add_trace(
                go.Violin(
                    y=vals,
                    name=method,
                    legendgroup=method,
                    showlegend=(col_idx == 1),
                    line_color=colors[method],
                    fillcolor=colors[method],
                    opacity=0.5,
                    box_visible=True,
                    meanline_visible=True,
                    points=False,
                ),
                row=1,
                col=col_idx,
            )
        fig.update_xaxes(title_text="Method", row=1, col=col_idx)
        if log_y:
            fig.update_yaxes(type="log", row=1, col=col_idx)
    y_title = f"{metric_title} (log scale)" if log_y else metric_title
    fig.update_yaxes(title_text=y_title, row=1, col=1)
    fig.update_layout(legend_title="Method", violinmode="group")
    return _apply_layout(fig)


def _build_grouped_run_box_figure(
    run_a_rows: Sequence[Dict[str, Any]],
    run_b_rows: Sequence[Dict[str, Any]],
    run_a_short_label: str,
    run_b_short_label: str,
    methods: Sequence[str],
    colors: Dict[str, str],
    metric: str,
    metric_title: str,
    log_y: bool = False,
    log_floor: float = MIN_POS_FLOAT,
) -> go.Figure:
    fig = go.Figure()
    for method_idx, method in enumerate(methods):
        base_color = colors.get(method, "#7f7f7f")
        run_specs = [
            (run_a_short_label, run_a_rows, _blend_with_white(base_color, 0.35), -0.2),
            (run_b_short_label, run_b_rows, _blend_with_black(base_color, 0.18), 0.2),
        ]
        for run_label, rows, color, x_delta in run_specs:
            vals = _filter_metric_values(rows, method, metric, log_axis=log_y, log_floor=log_floor)
            if vals.size == 0:
                continue
            x_pos = float(method_idx) + float(x_delta)
            fig.add_trace(
                go.Box(
                    x=[x_pos] * int(vals.size),
                    y=vals.tolist(),
                    name=f"{method} ({run_label})",
                    legendgroup=f"{method}_{run_label}",
                    width=0.32,
                    marker_color=color,
                    boxpoints=False,
                )
            )
    fig.update_layout(boxmode="overlay", legend_title="Method / benchmark")
    fig.update_xaxes(
        title_text="Method",
        tickmode="array",
        tickvals=[float(i) for i in range(len(methods))],
        ticktext=list(methods),
        range=[-0.6, max(float(len(methods)) - 0.4, 0.6)],
    )
    if log_y:
        fig.update_yaxes(type="log", title_text=f"{metric_title} (log scale)")
    else:
        fig.update_yaxes(title_text=metric_title)
    return _apply_layout(fig)


def _write_summary_csv(
    out_path: Path,
    run_rows: Sequence[Tuple[str, str, Sequence[Dict[str, Any]]]],
    methods: Sequence[str],
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "run_id",
        "run_label",
        "method",
        "total_cases",
        "success_cases",
        "success_rate_pct",
        "chamfer_median",
        "chamfer_p95",
        "hausdorff_median",
        "hausdorff_p95",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for run_id, run_label, rows in run_rows:
            for method in methods:
                all_rows = [r for r in rows if str(r.get("method", "")) == method]
                ok_rows = [r for r in all_rows if int(r.get("success", 0)) == 1]
                cham = np.asarray(
                    [float(r["chamfer_symmetric"]) for r in ok_rows if np.isfinite(float(r["chamfer_symmetric"]))],
                    dtype=np.float64,
                )
                haus = np.asarray(
                    [float(r["hausdorff"]) for r in ok_rows if np.isfinite(float(r["hausdorff"]))],
                    dtype=np.float64,
                )
                writer.writerow(
                    {
                        "run_id": run_id,
                        "run_label": run_label,
                        "method": method,
                        "total_cases": len(all_rows),
                        "success_cases": len(ok_rows),
                        "success_rate_pct": (100.0 * len(ok_rows) / len(all_rows)) if all_rows else float("nan"),
                        "chamfer_median": float(np.median(cham)) if cham.size else float("nan"),
                        "chamfer_p95": float(np.percentile(cham, 95)) if cham.size else float("nan"),
                        "hausdorff_median": float(np.median(haus)) if haus.size else float("nan"),
                        "hausdorff_p95": float(np.percentile(haus, 95)) if haus.size else float("nan"),
                    }
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate paper-focused plots from two benchmark runs.")
    parser.add_argument("--run-a", type=Path, default=DEFAULT_RUN_A, help="First run directory.")
    parser.add_argument("--run-b", type=Path, default=DEFAULT_RUN_B, help="Second run directory.")
    parser.add_argument("--run-a-label", type=str, default=None, help="Optional display label for run A.")
    parser.add_argument("--run-b-label", type=str, default=None, help="Optional display label for run B.")
    parser.add_argument("--run-a-short-label", type=str, default="L benchmark", help="Short legend label for run A.")
    parser.add_argument("--run-b-short-label", type=str, default="S benchmark", help="Short legend label for run B.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT, help="Output directory for paper plots.")
    parser.add_argument("--formats", type=str, default="pdf,html", help="Comma-separated formats, e.g. pdf,html.")
    parser.add_argument("--dpi", type=int, default=180, help="Export DPI for static image formats.")
    parser.add_argument(
        "--log-floor",
        type=float,
        default=DEFAULT_LOG_FLOOR,
        help="Lower clamp applied in log-based plots (values below are shown as this floor).",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_a = args.run_a.resolve()
    run_b = args.run_b.resolve()
    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    log_floor = max(float(args.log_floor), MIN_POS_FLOAT)

    formats = [s.strip().lower() for s in str(args.formats).split(",") if s.strip()]
    if not formats:
        formats = ["pdf", "html"]

    run_a_rows = _load_plot_rows(run_a)
    run_b_rows = _load_plot_rows(run_b)
    run_a_label = args.run_a_label or _regime_label_from_config(run_a)
    run_b_label = args.run_b_label or _regime_label_from_config(run_b)
    run_a_short_label = str(args.run_a_short_label)
    run_b_short_label = str(args.run_b_short_label)

    methods = sorted({str(r["method"]) for r in run_a_rows} | {str(r["method"]) for r in run_b_rows})
    palette = px.colors.qualitative.D3
    colors = {m: palette[i % len(palette)] for i, m in enumerate(methods)}

    created_files: List[str] = []
    warnings: List[str] = []

    plot_specs = [
        ("chamfer_ecdf", "chamfer_symmetric", "Chamfer distance (symmetric)", False, "ecdf"),
        ("chamfer_ecdf_logx", "chamfer_symmetric", "Chamfer distance (symmetric)", True, "ecdf"),
        ("chamfer_ecdf_grouped", "chamfer_symmetric", "Chamfer distance (symmetric)", True, "ecdf_grouped_chamfer"),
        ("chamfer_y_to_z_ecdf", "chamfer_y_to_z", "Chamfer distance (Y → Z)", False, "ecdf"),
        ("chamfer_y_to_z_ecdf_logx", "chamfer_y_to_z", "Chamfer distance (Y → Z)", True, "ecdf"),
        ("chamfer_z_to_y_ecdf", "chamfer_z_to_y", "Chamfer distance (Z → Y)", False, "ecdf"),
        ("chamfer_z_to_y_ecdf_logx", "chamfer_z_to_y", "Chamfer distance (Z → Y)", True, "ecdf"),
        ("chamfer_box_logy", "chamfer_symmetric", "Chamfer distance (symmetric)", True, "box"),
        ("chamfer_box_grouped_runs_logy", "chamfer_symmetric", "Chamfer distance (symmetric)", True, "grouped_box"),
        ("chamfer_violin", "chamfer_symmetric", "Chamfer distance (symmetric)", False, "violin"),
        ("chamfer_violin_logy", "chamfer_symmetric", "Chamfer distance (symmetric)", True, "violin"),
        ("hausdorff_ecdf", "hausdorff", "Hausdorff distance", False, "ecdf"),
        ("hausdorff_ecdf_logx", "hausdorff", "Hausdorff distance", True, "ecdf"),
        ("hausdorff_box_logy", "hausdorff", "Hausdorff distance", True, "box"),
        ("hausdorff_box_grouped_runs_logy", "hausdorff", "Hausdorff distance", True, "grouped_box"),
        ("hausdorff_violin", "hausdorff", "Hausdorff distance", False, "violin"),
        ("hausdorff_violin_logy", "hausdorff", "Hausdorff distance", True, "violin"),
    ]

    for name, metric, title, log_flag, kind in plot_specs:
        try:
            if kind == "ecdf":
                fig = _build_ecdf_figure(
                    run_a_rows,
                    run_b_rows,
                    run_a_label,
                    run_b_label,
                    methods,
                    metric,
                    title,
                    log_x=log_flag,
                    colors=colors,
                    log_floor=log_floor,
                )
            elif kind == "ecdf_grouped_chamfer":
                fig = _build_grouped_chamfer_ecdf_figure(
                    run_a_rows,
                    run_b_rows,
                    run_a_label,
                    run_b_label,
                    methods,
                    colors=colors,
                    log_x=log_flag,
                    log_floor=log_floor,
                )
            elif kind == "violin":
                fig = _build_violin_figure(
                    run_a_rows,
                    run_b_rows,
                    run_a_label,
                    run_b_label,
                    methods,
                    metric,
                    title,
                    colors=colors,
                    log_y=log_flag,
                    log_floor=log_floor,
                )
            elif kind == "grouped_box":
                fig = _build_grouped_run_box_figure(
                    run_a_rows,
                    run_b_rows,
                    run_a_short_label,
                    run_b_short_label,
                    methods,
                    colors,
                    metric,
                    title,
                    log_y=log_flag,
                    log_floor=log_floor,
                )
            else:
                fig = _build_log_box_figure(
                    run_a_rows,
                    run_b_rows,
                    run_a_label,
                    run_b_label,
                    methods,
                    metric,
                    title,
                    colors=colors,
                    log_floor=log_floor,
                )
            created_files.extend(_save_figure(fig, out_dir / name, formats=formats, dpi=int(args.dpi)))
        except Exception as exc:
            warnings.append(f"{name}: {exc}")

    _write_summary_csv(
        out_dir / "paper_summary_by_method.csv",
        [
            ("run_a", run_a_label, run_a_rows),
            ("run_b", run_b_label, run_b_rows),
        ],
        methods,
    )

    summary = {
        "run_a": str(run_a),
        "run_b": str(run_b),
        "run_a_label": run_a_label,
        "run_b_label": run_b_label,
        "run_a_short_label": run_a_short_label,
        "run_b_short_label": run_b_short_label,
        "output_dir": str(out_dir),
        "log_floor": log_floor,
        "methods": methods,
        "created_files": created_files,
        "warnings": warnings,
    }
    with (out_dir / "paper_plots_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))
    return 0 if not warnings else 2


if __name__ == "__main__":
    raise SystemExit(main())
