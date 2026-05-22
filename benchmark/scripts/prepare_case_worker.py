import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare a benchmark case and write JSON output.")
    parser.add_argument("--case-json", required=True, type=Path)
    parser.add_argument("--config-json", required=True, type=Path)
    parser.add_argument("--prep-dir", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    src_dir = Path(__file__).resolve().parents[1] / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    import benchmark_pipeline as bp

    case = json.loads(args.case_json.read_text(encoding="utf-8"))
    cfg = json.loads(args.config_json.read_text(encoding="utf-8"))
    result = bp.prepare_case_data(case, cfg, args.prep_dir)
    args.out_json.write_text(json.dumps(bp._to_serializable(result)), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
