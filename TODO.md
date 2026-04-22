# ephacRTools — To-Do List

---

## P0 · MultiSE Architecture (`multiSE` branch)

### ✅ M1 · Rename columns to snake_case throughout R package
**Scope:** `prepareSE.R`, `imagingTools.R`, `classificationTools.R`, `utilityTools.R`, `server.R`, `localImgBrowserModule.R`, tests, vignette

| Old | New |
|-----|-----|
| `Well` | `well_id` |
| `Plate_ID` | `plate_id` |
| `QC` | `qc` |
| `V_Clamp` | `v_clamp_mV` |
| `Sweep` (column) | `sweep` |
| `Row` / `Column` | `row` / `col` |

This is a **breaking change** — done once here on `multiSE` so it never needs revisiting. Coordinate with Python naming in `inst/python/`.

---

### ✅ M2 · Implement `MultiAssayExperiment` container
**New file:** `R/multiProtocol.R`

- `buildMAE(experiments, coldata)` — assembles a `MultiAssayExperiment` from a named list of SCE objects (one per protocol) + shared `colData`
- `sampleMap` auto-built from `well_id` + `plate_id` columns present in each experiment's `colData`
- `colData(mae)` holds shared well metadata (imaging, classification results, conditions)
- Each experiment retains its own `rowData` (protocol-specific sweep metadata)
- Add `MultiAssayExperiment` to `DESCRIPTION Imports`

---

### M3 · Update `mergeSEandImg()` and `mergeClassificationToSE()` for MAE
**Files:** `R/imagingTools.R`, `R/classificationTools.R`

Currently write to `colData(se)`. Add MAE-aware variants that write to `colData(mae)` instead, keyed by `well_id.plate_id`. Keep single-SE signatures unchanged for backwards compatibility.

---

### M4 · Update `reducedDim.Cellwise()` for multi-protocol
**File:** `R/utilityTools.R`

- Current: works on one SCE
- Add a wrapper that concatenates feature matrices across experiments (protocols) before PCA/UMAP, producing a joint embedding stored in `colData(mae)`
- Per-protocol embeddings remain available by calling the existing function on individual experiments

---

### M5 · Update `tinySEV` to handle MAE objects
**Files:** `R/ui.R`, `R/server.R`

- Detect whether the loaded object is SCE or MAE
- Add protocol selector (which experiment to display assay heatmaps / IV plots for)
- Joint embedding plots read from `colData(mae)` if available
- This is the largest UI change — tackle last after M1–M4 are stable

---

## P1 · Imaging SE (`multiSE` branch)

### I1 · `prepareImgSE()` — well-level imaging SCE from SQLite
**New file:** `R/imagingTools.R` (replaces/wraps current `prepareSingleImgDF` + `mergeSEandImg`)

Structure:
- **Rows** = channels (C1, C2, C3 — or user-labelled GFP, mCherry, BF)
- **Columns** = `well_id.plate_id` — same shape as electrophysiology SE, ready to join in MAE
- **Assays**: `mean_intensity`, `mean_area`, `total_area`, `norm_area`, `particle_count`, `total_intden`
  - aggregated from `allSelected` selection (hole-level ROI) in `PA_Measurement_Tables`
  - `norm_area = total_area / hole_area` computed from the `allSelected` `Selection_Area`
- **`rowData`**: `channel_name` (C1/C2/C3), `channel_label` (user-supplied mapping e.g. `c(C1="GFP", C2="mCherry")`)
- **`colData`**: `well_id`, `plate_id`, experimental metadata from `Particle_Analysis_Table`

```r
prepareImgSE(pathDB, channel_labels = c(C1 = "GFP", C2 = "mCherry", C3 = "BF"))
```

This replaces the current nested-DataFrame approach in `colData(se)`. Existing `mergeSEandImg()` becomes `addImgToMAE()` that attaches the imaging SCE as an experiment.

---

### I2 · `prepareParticleSE()` — particle-level SCE from SQLite
**New file:** `R/imagingTools.R`

Structure (single-cell analogy — particles = cells, morphological features = genes):
- **Rows** = morphological features: `Area`, `Mean`, `IntDen`, `Perim`, `Circ`, `AR`, `Round`, `Solidity` (~10 features)
- **Columns** = individual particles across all wells and channels
- **`colData`**: `well_id`, `plate_id`, `channel_name`, `particle_label`, `image_id`, `selection_type` (`hole` / `bf_image`)

```r
prepareParticleSE(pathDB, selection = c("hole", "bf_images", "all"))
```

Enables:
- `reducedDim.Cellwise()` on particle morphology → PCA/UMAP of particle shapes
- k-means / hierarchical clustering on morphology to find particle subtypes
- Aggregate cluster membership per well → classification (replaces manual `classificationTools.R` pipeline)
- `sampleMap` in MAE links each particle column back to its well in shared `colData`

---

### I3 · Update `classificationTools.R` to accept particle SCE
Currently takes a raw data frame. Add an SCE-aware path:
```r
classifyImgParticles(particle_se, ...)   # new: takes prepareParticleSE() output
classifyImgParticles(df, ...)            # old: still works for backwards compat
```
Classification results written to `colData(mae)` as flat columns.

---

### I4 · Update MAE assembly to include imaging experiments
**File:** `R/multiProtocol.R` (M2)

```r
mae <- buildMAE(
  ephys    = list(voltage_clamp = se_vc),
  imaging  = list(
    well     = prepareImgSE(db_path),
    particles = prepareParticleSE(db_path)
  )
)
```
`sampleMap` auto-generated: well-level experiments map 1:1, particle experiment maps n:1 via `well_id`.

---

## P2 · Current-Clamp / Python Pipeline (`multiSE` branch)

### ✅ C1 · Add `arrow` to `DESCRIPTION Imports`
One-liner. Required for reading parquet files without Python.

---

### ✅ C2 · Write `prepareCCSE()` from pre-computed sweep parquet
**New file:** `R/currentClampTools.R`

```r
prepareCCSE <- function(sweep_parquet_path, ap_parquet_path = NULL)
```

- Reads `sweep_table.parquet` (output of `ap_analysis_v2.py`) via `arrow::read_parquet()`
- Pivots into SCE: rows = sweeps, cols = `well_id.plate_id`
- Assays: `ap_count`, `ap_freq_hz`, `mean_ap_amplitude_mV`, `mean_isi_s`, `baseline_v_mV`, `steady_v_mV`, `mean_dvdt_max_mV_per_ms`, `mean_dvdt_min_mV_per_ms`
- `rowData`: `sweep`, `stim_amp_cmd`, `stim_start_s`, `stim_end_s`, `protocol`
- `colData`: `well_id`, `plate_id`, `date`
- If `ap_parquet_path` provided, attach per-spike data as nested `DataFrame` in `rowData`
- No Python dependency in this mode

---

### ✅ C3 · CSV → Parquet → AP analysis → SCE pipeline (redesigned as step functions)
**File:** `R/currentClampTools.R`

```r
prepareCCSE(path, mode = "raw", cfg = list())
```

- Calls `ap_analysis_v2.py::analyze_parquet_iterative()` via reticulate on the raw parquet
- Writes intermediate parquets to `tempdir()`, then calls C2 internally
- `reticulate` in `DESCRIPTION Suggests` (not `Imports`) — fail gracefully with a clear message if not installed
- Document required Python packages: `numpy`, `scipy`, `polars`, `pyarrow`

---

### ✅ C4 · Add `EfelConfig` R-side equivalent (`cc_config()`)
**File:** `R/currentClampTools.R`

```r
cc_config <- function(smooth_ms = 0.2, dvdt_thr = 25.0, refractory_ms = 2.0, ...)
```

Thin wrapper that builds the named list passed to the Python `EfelConfig` dataclass via reticulate. Lets users tune AP detection from R without touching Python.

---

### ✅ C5-extra · `plotCCTraces()` — raw trace overlay plot with stimulus inset and scalebars
**File:** `R/currentClampTools.R`
Supports: single well with inset, multi-well facets, viridis or fixed colour, scalebars, inset position control.

---

### ✅ C6-extra · `extractCCFeatures()` — per-well feature extraction from sweep/AP parquets
**File:** `R/currentClampTools.R`

Derives a one-row-per-well summary from a `prepareCCSE()` SCE (no raw traces needed):
- **Sweep parquet**: resting Vm, input resistance (MΩ), rheobase (pA), max firing rate (Hz), FI slope (Hz/pA), depolarisation block detection
- **AP parquet** (optional): first AP amplitude at rheobase, time to first spike (ms), AP amplitude tapering ratio, ISI adaptation ratio

---

### C5 · Write tests for `prepareCCSE()`
**File:** `tests/testthat/test-currentClamp.R`

Use a small fixture parquet (generated from the existing notebook data) bundled under `inst/testdata/`.

---

## P3 · Existing Bug Fixes (carry over from main, apply on `multiSE` branch)

### #6 · `addImgPaths()` — hardcoded path
**File:** `R/imagingTools.R` — replace hardcoded `pathImg` with the `parent_folder` parameter.

### #9 · `prepareSE()` — non-unique colnames across plates
**File:** `R/prepareSE.R` — keep colnames as `interaction(well_id, plate_id)`.

### #10 · Remove debug `print()`/`cat()` calls
Scattered across `prepareSE.R`, `imagingTools.R`, `utilityTools.R`.

### #11 · Remove `rm(ls())` / `gc()` from `prepareDF()`
Dangerous inside a function — remove both lines.

### #3 · Remove `print(c)` from `ag()`
**File:** `R/misc.R`.

---

## P4 · Infrastructure

### #12 · Fix `DESCRIPTION`
- Fill in `Description:` field
- Set `License:` (suggest MIT)
- Remove circular `Remotes: dcolam/ephacRTools`
- Add missing `Imports`: `jsonlite`, `rlang`, `arrow`

### #17 · Fix vignette
Replace `devtools::load_all()` and relative paths with `library(ephacRTools)` + bundled datasets.

### #18 · Fix `tryCatch` not catching `req()` in `server.R`
Add `shiny.silent.error` condition or guard with `isTruthy()`.

### #19 · Remove dead `.extract_plate_from_path()` from `server.R`

---

## P5 · Future / DataControl-Agnostic Pipeline

> Not for now — revisit once multiSE and current-clamp work is stable.

### F1 · DataControl-agnostic VC pipeline from raw `.json` + `.dat` files
**Feasibility: confirmed and benchmarked.** Format fully decoded from `inst/test-data/18T39412/` using `inst/test-data/python_read_trace_from_json.python`.

**Binary format (corrected):**
- **int16 little-endian**, NOT float32. Each value × per-well `I2DScale` factor = Amperes.
- Layout: `[sweep % SweepsPerFile][col][row][NofSamples × LeakData]` — pure, no header.
- `LeakData=2`: raw trace (samples 0..N) + **leak-subtracted trace** (samples N+1..2N) both stored — leak correction is free from the instrument.
- File size confirmed: `10 × 24 × 16 × (3100 × 2 × 2 bytes) = 47,616,000 bytes ✓`

**JSON keys:**
- `TraceHeader.TimeScalingIV.I2DScale` — 384-element per-well scale factors (IV protocols)
- `TraceHeader.TimeScaling.I2DScale` — same for regular protocols
- `TraceHeader.TimeScalingIV.Stimulus[sweep]` — per-sweep voltage command (IV)
- `ExperimentConditions.OAFunctions` — cursor list: `{CursorName, TimeStart_ms, TimeEnd_ms, EvalMethod_Enum}` (0=Min, 1=Max, 3=Integral, 4=Mean, 11=Erev)
- `QCData.RSeal/Capacitance/Rseries` — pre-computed `[sweep][col][row]` arrays

**Performance benchmark (pure Python, no numpy):**
- 11 sweeps × 384 wells × 2 cursors: **2.4 seconds**
- With numpy (vectorised file load): estimated **~50ms**
- Full recording day (5 protocols, ~400 sweeps): **<2 seconds** with numpy

**Custom OA functions:** trivially supported. JSON cursors reproduce DataControl built-ins automatically; `extra_cursors` parameter adds arbitrary user-defined analysis:
```r
prepareRawSE(folder,
  extra_cursors = list(
    list(name="tau_inact", t_start=80, t_end=200,
         fn = function(seg, t) fit_exp_tau(seg, t)),
    list(name="peak_vs_baseline", t_start=60, t_end=80,
         fn = function(seg, baseline) min(seg) - mean(baseline))
  )
)
```
Each cursor becomes a named assay in the output SCE — completely protocol-agnostic.

**Implementation plan (`inst/python/read_raw_syncropatch.py` + `R/prepareRawSE.R`):**
1. `parse_json_metadata()` — OAFunctions, I2DScale, sweep→compound map, QC arrays
2. `read_all_wells()` — load full `.dat` file as numpy int16 array, reshape `[sweeps, cols, rows, 2*samples]`, multiply by I2DScale
3. `apply_cursors()` — vectorised slice + reduce per cursor over `[sweeps, cols, rows]`
4. `assemble_sce()` — same output shape as `prepareSE()`

Reference script: `inst/test-data/python_read_trace_from_json.python` (Nanion-provided, single-well reader — basis for the all-wells vectorised version).

### F2 · Unified `prepareProtocol()` entry point
Once both Excel-based VC and parquet-based CC pipelines exist, a single dispatcher:
```r
prepareProtocol(path, protocol = c("voltage_clamp", "current_clamp"), ...)
```
that routes to `prepareSE()` or `prepareCCSE()` and returns an SCE ready to drop into `buildMAE()`.
