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
  "save_cut_meshes": false,
  "save_expected_result": false
}
```

- `save_cut_meshes=false` (recommended for large runs): do not persist intermediate cut meshes (`A/B/X/C/D/Z`, cutters, cut steps).
- `save_cut_meshes=true`: keep all intermediate meshes for visual debugging.
- `save_expected_result=true`: save expected output mesh `Z.obj` per case even when `save_cut_meshes=false`.
- Expected/debug meshes are optionally written once per case into `<output-dir>\expected_results\case_xxxxx\` (only when corresponding debug flags are enabled).
- In multi-method runs, expensive preparation (`X`, `C`, `D`, `Z`) is computed once per case and reused for all methods.

## Time and memory limits

`default_config.json` contains execution limits for both case preparation and tested methods:

```json
"limits": {
  "preparation_timeout_sec": 7200,
  "preparation_memory_limit_mb": 32768,
  "method_timeout_sec": 7200,
  "method_memory_limit_mb": 32768,
  "memory_check_interval_sec": 0.01
}
```

- `preparation_timeout_sec`: max time for case preparation (`X/C/D/Z` construction).
- `preparation_memory_limit_mb`: max RSS memory for preparation worker process.
- `method_timeout_sec`: max time per method execution.
- `method_memory_limit_mb`: max RSS memory per method process.
- `memory_check_interval_sec`: process polling interval (same style as your `BooleansTest.py` monitor loop).

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

## Plot generation

`default_config.json` now supports:

```json
"plots": {
  "enabled": false,
  "output_subdir": "plots",
  "dpi": 150,
  "formats": ["pdf"],
  "complexity_bins": 80,
  "export_timeout_sec": 30.0,
  "fallback_html_on_export_failure": true
}
```

When `plots.enabled=true`, the run automatically creates:
- `hausdorff_by_method`
- `chamfer_by_method`
- `runtime_vs_input_size` (elapsed time vs combined input size in MB)
- `memory_vs_input_size` (peak memory vs combined input size in MB)

Plot implementation is based on **Plotly** (matching the style used in `EvaluationScripts\scatter_plots.py` and `success_graph.py`).
If static export stalls/fails in your environment, the pipeline times out static image export and writes `.html` plot files as fallback.

Plot-only mode (no benchmark execution):

```powershell
python benchmark\scripts\run_pipeline.py `
  --config benchmark\config\default_config.json `
  --output-dir C:\path\to\existing_run `
  --plots-only
```

This reads existing method outputs from `--output-dir` and writes plots into `<output-dir>\plots\`.

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
- Method execution metadata includes `runtime_sec`, `peak_rss_mb`, and combined input size (`input_size_mb`) computed from method inputs (`C + D`) in MB, used by complexity plots.
