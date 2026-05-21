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
  "save_output_meshes": true
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

## Debug output control

`benchmark\config\default_config.json` contains:

```json
"debug": {
  "save_cut_meshes": false
}
```

- `save_cut_meshes=false` (recommended for large runs): do not persist intermediate cut meshes (`A/B/X/C/D/Z`, cutters, cut steps).
- `save_cut_meshes=true`: keep all intermediate meshes for visual debugging.

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

Dataset mode now:
- computes and caches dataset mesh stats in `<dataset>\mesh_stats.json` and `<dataset>\mesh_stats.csv`
- samples cases only from meshes with `is_closed=true` and `is_manifold=true`
- for STL meshes, topology checks use welded vertices (duplicate STL vertices are merged for closed/manifold detection)
- if fewer valid unique input pairs exist than requested via `--pairs`, it runs all available pairs
- writes selected pairs to `<method output>\input_mesh_pairs.json/csv`

Use `--update-stats` to force mesh stats recomputation:

```powershell
python benchmark\scripts\run_pipeline.py `
  --methods-dir benchmark\config\methods `
  --dataset-dir C:\path\to\dataset `
  --pairs 64 `
  --update-stats `
  --output-dir C:\path\to\output
```

## Key outputs per case

- `case_result.json` (full metadata, cut logs, metrics)
- Intermediate cut meshes are saved only when `debug.save_cut_meshes=true`.
- Result meshes are saved in `<method output_path>\meshes\` when `save_output_meshes=true`.

