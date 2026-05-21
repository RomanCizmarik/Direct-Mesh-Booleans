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
  --config benchmark\config\default_config.json `
  --methods-dir benchmark\config\methods `
  --input-a C:\path\to\A.obj `
  --input-b C:\path\to\B.obj `
  --op union `
  --output-dir C:\path\to\run_root `
  --case-name my_case
```

## Per-method configuration

Method runner settings are defined **per binary** via JSON files in `benchmark\config\methods\`.
The pipeline loads all method configs from that folder (or one file via `--method-config`).

Config template:

```json
{
  "exe_path": "C:\\path\\to\\method.exe",
  "operation_map": {
    "union": "...",
    "diff": "...",
    "int": "..."
  },
  "output_path": "method_output_dir",
  "save_output_meshes": true,
  "allow_empty_output_mesh": false
}
```

- `exe_path`: tested binary path
- `operation_map`: operation argument mapping
  - internal op `union` maps to `operation_map["union"]`
  - internal op `difference` maps to `operation_map["diff"]`
  - internal op `intersection` maps to `operation_map["int"]`
- `output_path`: where method-specific benchmark results are written
  - absolute path => used directly
  - relative path => relative to `--output-dir`
- `save_output_meshes`: if `true`, saves produced `Y` meshes to `<output_path>\meshes\`
- `allow_empty_output_mesh` (optional, default `false`): if `false`, header-only/empty meshes are treated as method failure

Provided configs:
- `direct_mesh_booleans.json`
- `direct_mesh_booleans_extension.json`
- `geogram.json`
- `volume_mesher.json`
- `libigl.json`

Example (single method):

```powershell
python benchmark\scripts\run_pipeline.py `
  --config benchmark\config\default_config.json `
  --method-config benchmark\config\methods\direct_mesh_booleans.json `
  --input-a C:\skola\PhD\Direct-Mesh-Booleans\data\sphere1.obj `
  --input-b C:\skola\PhD\Direct-Mesh-Booleans\data\sphere2.obj `
  --op union `
  --output-dir benchmark\artifacts\sphere_test_dmb_cfg `
  --case-name sphere_pair_test
```

Example (all methods in `config\methods`):

```powershell
python benchmark\scripts\run_pipeline.py `
  --config benchmark\config\default_config.json `
  --methods-dir benchmark\config\methods `
  --input-a C:\skola\PhD\Direct-Mesh-Booleans\data\sphere1.obj `
  --input-b C:\skola\PhD\Direct-Mesh-Booleans\data\sphere2.obj `
  --op union `
  --output-dir benchmark\artifacts\sphere_test_methods_template `
  --case-name sphere_pair_test
```

## Run random dataset benchmark

```powershell
python benchmark\scripts\run_pipeline.py `
  --methods-dir benchmark\config\methods `
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

