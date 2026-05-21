from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path


REPO_BENCHMARK_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_BENCHMARK_DIR / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from benchmark_pipeline import (  # noqa: E402
    _normalize_op,
    load_config,
    load_method_configs,
    run_dataset_benchmark,
    run_manual_case,
    save_global_outputs,
)


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
    parser.add_argument("--input-a", type=Path, default=None, help="Manual case: first mesh path.")
    parser.add_argument("--input-b", type=Path, default=None, help="Manual case: second mesh path.")
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

    root_out_dir = args.output_dir.resolve()
    root_out_dir.mkdir(parents=True, exist_ok=True)

    per_method_overview = []

    for method_cfg in methods:
        method_name = str(method_cfg.get("method_name", "method"))
        method_output_path = Path(str(method_cfg["output_path"]))
        if not method_output_path.is_absolute():
            method_output_path = root_out_dir / method_output_path
        method_output_path.mkdir(parents=True, exist_ok=True)

        run_cfg = copy.deepcopy(config)
        run_cfg["method_under_test"] = method_cfg

        if args.input_a is not None or args.input_b is not None:
            if args.input_a is None or args.input_b is None:
                raise ValueError("Both --input-a and --input-b must be provided for manual case mode.")
            operation = _normalize_op(args.op)
            result = run_manual_case(
                input_a=args.input_a.resolve(),
                input_b=args.input_b.resolve(),
                operation=operation,
                case_name=args.case_name,
                config=run_cfg,
                output_root=method_output_path,
            )
            save_global_outputs(method_output_path, run_cfg, [result])
            per_method_overview.append(
                {
                    "method": method_name,
                    "status": result.get("status"),
                    "case_id": result.get("case_id"),
                    "output_dir": str(method_output_path),
                }
            )
            continue

        if args.dataset_dir is None:
            raise ValueError("Dataset mode requires --dataset-dir unless manual inputs are provided.")

        results = run_dataset_benchmark(args.dataset_dir.resolve(), run_cfg, method_output_path)
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
    print(json.dumps(per_method_overview, indent=2))
    has_failures = any(item.get("status") == "error" or item.get("failed", 0) > 0 for item in per_method_overview)
    return 2 if has_failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

