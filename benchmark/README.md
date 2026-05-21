# Broken-input Boolean benchmark (Python)

This folder contains a Python benchmark pipeline for:
1. closed-input reference Boolean `X` via libigl,
2. sphere-based breakage of inputs to produce `C` and `D`,
3. expected output `Z` by source-aware cuts on `X`,
4. method-under-test output `Y` on `(C,D)`,
5. metrics and per-case artifacts.

## Install

```powershell
python -m pip install -r benchmark\requirements.txt
```

## Run a manual pair

```powershell
python benchmark\scripts\run_pipeline.py `
  --input-a C:\path\to\A.obj `
  --input-b C:\path\to\B.obj `
  --op union `
  --output-dir C:\path\to\output `
  --case-name my_case
```

## Run random dataset benchmark

```powershell
python benchmark\scripts\run_pipeline.py `
  --dataset-dir C:\path\to\dataset `
  --pairs 64 `
  --output-dir C:\path\to\output
```

## Key outputs per case

- `A.obj`, `B.obj` (original inputs)
- `X.obj` + `X_birth_indices.npy` (closed-input reference and provenance)
- `cutters\A\sphere_*.obj`, `cutters\B\sphere_*.obj` (all generated cutting spheres)
- `C.obj`, `D.obj` (broken inputs)
- `Z.obj` (expected broken-output reference)
- `Y.obj` (method-under-test output)
- `case_result.json` (full metadata, cut logs, metrics)

