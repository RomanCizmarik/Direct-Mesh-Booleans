from __future__ import annotations

import argparse
import copy
import csv
import json
import shutil
import sys
from pathlib import Path

import numpy as np


REPO_BENCHMARK_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_BENCHMARK_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from benchmark_pipeline import (  # noqa: E402
    _normalize_op,
    get_dataset_mesh_stats,
    load_config,
    load_mesh,
    load_method_configs,
    prepare_case_data_with_limits,
    run_case_with_prepared,
    save_global_outputs,
)
from plots import generate_standard_plots  # noqa: E402


def _write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _write_csv(path: Path, rows: list[dict]) -> None:
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
            writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Boolean benchmark pipeline for broken mesh inputs.")
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_BENCHMARK_DIR / "config" / "default_config.json",
        help="Path to JSON config.",
    )
    parser.add_argument(
        "--method-config",
        type=Path,
        default=None,
        help="Path to one method JSON config. If omitted, all configs in --methods-dir are used.",
    )
    parser.add_argument(
        "--methods-dir",
        type=Path,
        default=REPO_BENCHMARK_DIR / "config" / "methods",
        help="Directory with method JSON configs (used when --method-config is omitted).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_BENCHMARK_DIR / "artifacts" / "run",
        help="Output directory for manifests, meshes, and results.",
    )
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        default=None,
        help="Dataset directory (used when --input-a/--input-b are not provided).",
    )
    parser.add_argument("--pairs", type=int, default=None, help="Override number of sampled dataset pairs.")
    parser.add_argument("--seed", type=int, default=None, help="Override global seed.")
    parser.add_argument(
        "--update-stats",
        action="store_true",
        help="Recompute dataset mesh stats cache (mesh_stats.json/csv) before running.",
    )
    parser.add_argument("--input-a", type=Path, default=None, help="Manual case: first mesh path.")
    parser.add_argument("--input-b", type=Path, default=None, help="Manual case: second mesh path.")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Generate plots only from --output-dir using existing results (no benchmark execution).",
    )
    parser.add_argument(
        "--stats-only",
        action="store_true",
        help="Run only dataset mesh stats stage and write mesh_stats.json/csv.",
    )
    parser.add_argument(
        "--op",
        type=str,
        default="union",
        help="Manual case operation: union|intersection|difference (or aliases).",
    )
    parser.add_argument("--case-name", type=str, default="manual_case", help="Manual case identifier.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config(args.config)
    root_out_dir = args.output_dir.resolve()
    root_out_dir.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        plots_cfg = copy.deepcopy(config.get("plots", {}))
        plots_cfg["enabled"] = True
        plot_summary = generate_standard_plots(root_out_dir, plots_cfg)
        print(json.dumps({"plots": plot_summary}, indent=2))
        return 0
    if args.stats_only:
        if args.dataset_dir is None:
            raise ValueError("--stats-only requires --dataset-dir.")
        dataset_dir = args.dataset_dir.resolve()
        mesh_records = get_dataset_mesh_stats(dataset_dir, update_stats=True)
        print(
            json.dumps(
                {
                    "dataset_dir": str(dataset_dir),
                    "mesh_count": len(mesh_records),
                    "mesh_stats_json": str(dataset_dir / "mesh_stats.json"),
                    "mesh_stats_csv": str(dataset_dir / "mesh_stats.csv"),
                },
                indent=2,
            )
        )
        return 0

    if args.pairs is not None:
        config["num_pairs"] = int(args.pairs)
    if args.seed is not None:
        config["seed"] = int(args.seed)

    methods: list[dict]
    if args.method_config is not None:
        methods = load_method_configs(args.method_config.parent.resolve())
        methods = [m for m in methods if Path(m["config_path"]).resolve() == args.method_config.resolve()]
        if not methods:
            raise ValueError(f"Method config not found: {args.method_config}")
    else:
        methods = load_method_configs(args.methods_dir.resolve())

    expected_results_dir = root_out_dir / "expected_results"
    global_debug_cfg = config.get("debug", {})
    debug_requested = bool(global_debug_cfg.get("save_cut_meshes", False) or global_debug_cfg.get("save_expected_result", False))

    method_runs: list[dict] = []
    for method_cfg in methods:
        method_name = str(method_cfg.get("method_name", "method"))
        method_output_path = Path(str(method_cfg["output_path"]))
        if not method_output_path.is_absolute():
            method_output_path = root_out_dir / method_output_path
        method_output_path.mkdir(parents=True, exist_ok=True)

        run_cfg = copy.deepcopy(config)
        run_cfg["method_under_test"] = method_cfg
        run_cfg.setdefault("debug", {})
        if debug_requested:
            run_cfg["debug"]["expected_results_dir"] = str(expected_results_dir)
        else:
            run_cfg["debug"].pop("expected_results_dir", None)

        method_runs.append(
            {
                "method_name": method_name,
                "method_cfg": method_cfg,
                "output_path": method_output_path,
                "run_cfg": run_cfg,
                "results": [],
            }
        )

    manual_mode = args.input_a is not None or args.input_b is not None
    manual_cases: list[dict] = []
    candidate_records: list[dict] = []
    target_case_count = 0
    dataset_pair_indices: list[tuple[int, int]] = []
    dataset_rng = np.random.default_rng(int(config["seed"]))
    normalized_ops = [_normalize_op(op) for op in list(config["operations"])]
    if not normalized_ops:
        raise ValueError("At least one operation must be provided.")

    if manual_mode:
        if args.input_a is None or args.input_b is None:
            raise ValueError("Both --input-a and --input-b must be provided for manual case mode.")
        operation = _normalize_op(args.op)
        manual_cases = [
            {
                "case_id": args.case_name,
                "input_a": str(args.input_a.resolve()),
                "input_b": str(args.input_b.resolve()),
                "operation": operation,
                "seed_spheres_a": int(dataset_rng.integers(0, 2**31 - 1)),
                "seed_spheres_b": int(dataset_rng.integers(0, 2**31 - 1)),
            }
        ]
        target_case_count = len(manual_cases)
    else:
        if args.dataset_dir is None:
            raise ValueError("Dataset mode requires --dataset-dir unless manual inputs are provided.")
        dataset_dir = args.dataset_dir.resolve()
        mesh_records = get_dataset_mesh_stats(dataset_dir, update_stats=bool(args.update_stats))
        dataset_manifest_dir = root_out_dir / "manifests"
        _write_json(dataset_manifest_dir / "mesh_manifest.json", mesh_records)
        _write_csv(dataset_manifest_dir / "mesh_manifest.csv", mesh_records)
        candidate_records = [
            r
            for r in mesh_records
            if r.get("status") == "ok"
            and bool(r.get("is_closed", False))
            and bool(r.get("is_manifold", False))
            and bool(r.get("is_winding_consistent", False))
            and not bool(r.get("inside_out_suspect", True))
        ]
        dataset_pair_indices = [(i, j) for i in range(len(candidate_records)) for j in range(i + 1, len(candidate_records))]
        dataset_rng.shuffle(dataset_pair_indices)
        target_case_count = min(int(config["num_pairs"]), len(dataset_pair_indices))

    prep_cfg = copy.deepcopy(config)
    prep_cfg.setdefault("debug", {})
    if debug_requested:
        prep_cfg["debug"]["expected_results_dir"] = str(expected_results_dir)
    else:
        prep_cfg["debug"].pop("expected_results_dir", None)

    prep_root_dir = root_out_dir / "_prepared_cases"
    cases: list[dict] = []
    input_pairs: list[dict] = []

    if manual_mode:
        for case in manual_cases:
            prep_case_dir = prep_root_dir / str(case["case_id"])
            if prep_case_dir.exists():
                shutil.rmtree(prep_case_dir, ignore_errors=True)
            prepared = prepare_case_data_with_limits(case, prep_cfg, prep_case_dir)
            if prepared.get("status") == "ok":
                v_z, f_z = load_mesh(Path(str(prepared["z_path"])))
                prepared["v_z"] = v_z
                prepared["f_z"] = f_z
            for run in method_runs:
                result = run_case_with_prepared(case, run["run_cfg"], run["output_path"], prepared)
                run["results"].append(result)
            if prep_case_dir.exists():
                shutil.rmtree(prep_case_dir, ignore_errors=True)
            cases.append(case)
            input_pairs.append(
                {
                    "case_id": case["case_id"],
                    "input_a": case["input_a"],
                    "input_b": case["input_b"],
                }
            )
    else:
        pair_cursor = 0
        while len(cases) < target_case_count and pair_cursor < len(dataset_pair_indices):
            ia, ib = dataset_pair_indices[pair_cursor]
            pair_cursor += 1
            case = {
                "case_id": f"case_{len(cases):05d}",
                "input_a": candidate_records[ia]["path"],
                "input_b": candidate_records[ib]["path"],
                "operation": normalized_ops[int(dataset_rng.integers(0, len(normalized_ops)))],
                "seed_spheres_a": int(dataset_rng.integers(0, 2**31 - 1)),
                "seed_spheres_b": int(dataset_rng.integers(0, 2**31 - 1)),
            }

            prep_case_dir = prep_root_dir / str(case["case_id"])
            if prep_case_dir.exists():
                shutil.rmtree(prep_case_dir, ignore_errors=True)
            prepared = prepare_case_data_with_limits(case, prep_cfg, prep_case_dir)

            valid_expected = False
            if prepared.get("status") == "ok":
                try:
                    v_z, f_z = load_mesh(Path(str(prepared["z_path"])))
                    if v_z.shape[0] > 0 and f_z.shape[0] > 0:
                        prepared["v_z"] = v_z
                        prepared["f_z"] = f_z
                        valid_expected = True
                except Exception:
                    valid_expected = False

            if not valid_expected:
                if prep_case_dir.exists():
                    shutil.rmtree(prep_case_dir, ignore_errors=True)
                continue

            for run in method_runs:
                result = run_case_with_prepared(case, run["run_cfg"], run["output_path"], prepared)
                run["results"].append(result)
            if prep_case_dir.exists():
                shutil.rmtree(prep_case_dir, ignore_errors=True)

            cases.append(case)
            input_pairs.append(
                {
                    "case_id": case["case_id"],
                    "input_a": case["input_a"],
                    "input_b": case["input_b"],
                }
            )

    _write_json(root_out_dir / "cases_manifest.json", cases)
    _write_csv(root_out_dir / "cases_manifest.csv", cases)
    _write_json(root_out_dir / "input_mesh_pairs.json", input_pairs)
    _write_csv(root_out_dir / "input_mesh_pairs.csv", input_pairs)

    per_method_overview = []
    for run in method_runs:
        method_name = str(run["method_name"])
        method_output_path = Path(run["output_path"])
        run_cfg = run["run_cfg"]
        results = run["results"]
        _write_json(method_output_path / "cases_manifest.json", cases)
        _write_csv(method_output_path / "cases_manifest.csv", cases)
        _write_json(method_output_path / "input_mesh_pairs.json", input_pairs)
        _write_csv(method_output_path / "input_mesh_pairs.csv", input_pairs)
        save_global_outputs(method_output_path, run_cfg, results)
        failed = sum(1 for r in results if r.get("status") != "ok")
        per_method_overview.append(
            {
                "method": method_name,
                "cases": len(results),
                "failed": failed,
                "output_dir": str(method_output_path),
            }
        )

    summary_path = root_out_dir / "methods_run_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(per_method_overview, handle, indent=2)

    plots_summary = None
    if bool(config.get("plots", {}).get("enabled", False)):
        plots_summary = generate_standard_plots(root_out_dir, config.get("plots", {}))
    print(json.dumps(per_method_overview, indent=2))
    if plots_summary is not None:
        print(json.dumps({"plots": plots_summary}, indent=2))
    has_failures = any(item.get("status") == "error" or item.get("failed", 0) > 0 for item in per_method_overview)
    return 2 if has_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
