from __future__ import annotations

import argparse
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
        "--method",
        type=str,
        default=None,
        choices=["libigl", "external_command"],
        help="Override method under test.",
    )
    parser.add_argument(
        "--method-command",
        type=str,
        default=None,
        help="External command template when method is external_command.",
    )
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
    if args.method is not None:
        config["method_under_test"]["name"] = args.method
    if args.method_command is not None:
        config["method_under_test"]["external_command_template"] = args.method_command

    out_dir = args.output_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.input_a is not None or args.input_b is not None:
        if args.input_a is None or args.input_b is None:
            raise ValueError("Both --input-a and --input-b must be provided for manual case mode.")
        operation = _normalize_op(args.op)
        result = run_manual_case(
            input_a=args.input_a.resolve(),
            input_b=args.input_b.resolve(),
            operation=operation,
            case_name=args.case_name,
            config=config,
            output_root=out_dir,
        )
        save_global_outputs(out_dir, config, [result])
        print(json.dumps({"status": result.get("status"), "case_id": result.get("case_id")}, indent=2))
        return 0 if result.get("status") == "ok" else 2

    if args.dataset_dir is None:
        raise ValueError("Dataset mode requires --dataset-dir unless manual inputs are provided.")

    results = run_dataset_benchmark(args.dataset_dir.resolve(), config, out_dir)
    save_global_outputs(out_dir, config, results)
    failed = sum(1 for r in results if r.get("status") != "ok")
    print(json.dumps({"cases": len(results), "failed": failed}, indent=2))
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

