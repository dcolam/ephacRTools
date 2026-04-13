# ephacRTools — Functions & Pipelines Reference

Version: **0.3.0** | Generated: 2026-04-01

---

## Table of Contents

1. [Package Overview](#1-package-overview)
2. [Data Ingestion Pipeline](#2-data-ingestion-pipeline)
3. [Imaging Pipeline](#3-imaging-pipeline)
4. [Classification Pipeline](#4-classification-pipeline)
5. [Clustering Pipeline](#5-clustering-pipeline)
6. [Utility & Analysis Functions](#6-utility--analysis-functions)
7. [Local Image Browser Module](#7-local-image-browser-module)
8. [Shiny App — tinySEV](#8-shiny-app--tinysev)
9. [Misc / Internal Helpers](#9-misc--internal-helpers)
10. [SCE Object Structure](#10-sce-object-structure)
11. [Known Issues](#11-known-issues)

---

## 1. Package Overview

`ephacRTools` processes **High-Throughput Automated Patch Clamp (HT-APC)** data from Syncropatch/Nanion instruments. The package converts DataControl Excel exports into `SingleCellExperiment` objects, integrates fluorescence imaging from Cluster Analysis SQLite databases, classifies cells by fluorescence phenotype, and exposes an interactive Shiny dashboard (`tinySEV`).

**Full data flow:**

```
DataControl Excel (.xlsx)
  └─► prepareDF / prepareMultipleDFs ──► prepareSE ──► SingleCellExperiment (SCE)

Cluster Analysis SQLite (.db)
  └─► previewImgDB (optional preview)
  └─► prepareSingleImgDF / prepareImgDF ──► df_cleaned ──► mergeSEandImg ──► SCE + nested colData

JPG thumbnail folder
  └─► addThumbnailPaths ──► SCE colData$thumbnails
      └─► plotThumbnails

Particle imaging data
  └─► filterParticles ──► aggregateParticles ──► scoreParticles ──► classifyWells
      (all-in-one: classifyImgParticles)
  └─► mergeClassificationToSE ──► SCE colData$Classification

SCE
  ├─► colAG / get_metric / fit_boltzmann_se   (electrophysiology metrics)
  ├─► reducedDim.Cellwise / plotDimRed        (dimensionality reduction)
  ├─► subsetSE / checkAssay                   (filtering / validation)
  └─► tinySEV()                               (interactive Shiny dashboard)
```

---

## 2. Data Ingestion Pipeline

**File:** `R/prepareSE.R`

### `prepareDF(pathToDF)`

Reads a single DataControl Excel file and returns a **long-format** `data.frame`.

| Argument | Type | Description |
|---|---|---|
| `pathToDF` | character | Path to a `.xlsx` file exported from DataControl |

**Returns:** `data.frame` with columns `Well`, `QC`, `Plate_ID`, one column per Online Analysis measurement, `Sweep`, and `V_Clamp`.

**What it does:**
- Reads the sheet whose name contains `"Export"`
- Standardises well names and plate barcode (`Nanion Chip Barcode` → `Plate_ID`)
- Extracts per-sweep voltage from `"Sweep Voltage"` rows when present
- Pivots sweep-wise columns into long format (one row per well × sweep)
- Auto-detects column types via `hablar::retype()`

**Known issues:** Contains `cat()` / `print()` debug statements and a `rm(list = setdiff(ls(), "new.df"))` call that should be removed (see TODO #10, #11).

---

### `prepareMultipleDFs(pathList)`

Calls `prepareDF()` on each path and `dplyr::bind_rows()` the results.

| Argument | Type | Description |
|---|---|---|
| `pathList` | character vector | Paths to one or more `.xlsx` files |

**Returns:** Single long-format `data.frame` with an additional `column_label` column (source filename).

---

### `prepareSE(pathDF, conditionColumns = c("Compound"))`

Main entry point. Reads one or more Excel files and builds a `SingleCellExperiment`.

| Argument | Type | Default | Description |
|---|---|---|---|
| `pathDF` | character (vector) | — | Path(s) to DataControl `.xlsx` file(s) |
| `conditionColumns` | character vector | `c("Compound")` | Columns describing experimental conditions |

**Returns:** `SingleCellExperiment` with:
- **assays** — one named matrix per numeric measurement (rows = sweeps, cols = wells)
- **rowData** — `Sweep`; per-sweep voltage / compound info (collapsed if uniform, nested DataFrame if variable across wells)
- **colData** — `Well`, `QC`, `Plate_ID`, `Row`, `Column`
- `colnames(se)` — `interaction(Well, Plate_ID)` format (e.g. `"A01.Plate1"`)

**Known issue:** Bug #9 — colnames are correctly set as `interaction(Well, Plate_ID)` but an earlier code path used `cd$Well` alone, risking non-unique names across plates. Verify this is resolved for multi-plate imports.

---

## 3. Imaging Pipeline

**File:** `R/imagingTools.R`

### `previewImgDB(paths, analysis = "pa")`

Scans one or more SQLite databases **before** import to discover available columns, channels, selection values, and Image_ID ranks. Used by `tinySEV` to populate the Import Data UI.

| Argument | Default | Description |
|---|---|---|
| `paths` | — | Character vector of `.db` file paths |
| `analysis` | `"pa"` | `"pa"` (Particle Analysis) or `"coloc"` (Colocalization) |

**Returns:** Named list with `pa_cols`, `meas_cols`, `selections`, `channels`, `n_rows`, `n_wells`, `image_ranks`, and per-DB summaries.

---

### `prepareSingleImgDF(pathDB, ...)`

Imports particle-level data from a single Cluster Analysis SQLite database.

| Key arguments | Description |
|---|---|
| `pathDB` | Path to `.db` file |
| `analysis` | `"pa"` or `"coloc"` |
| `selections` | Which Selection values to keep (e.g. `c("In-Focus")`) |
| `image_rank` | Which Image_ID rank to keep per well (1 = lowest, typically fluorescence) |
| `meas_cols` | Measurement columns to import |

**Returns:** `data.frame` with `Well`, `Plate_ID`, `Channel_Name`, `Selection`, `Area`, `Mean`, `IntDen`, `normArea`, plus requested measurement columns.

---

### `prepareImgDF(pathDB, ...)`

Wrapper for multi-database import; calls `prepareSingleImgDF()` on each file and `bind_rows()` the results.

---

### `df_cleaned(df, channels = c("Green", "Red", "ROMK"))`

Labels ROIs as `"Hole_ROI"` (smallest area per selection group) or `"background_ROI"` (largest area), adds a `CorrSel` column.

| Argument | Default | Description |
|---|---|---|
| `df` | — | Output of `prepareImgDF()` |
| `channels` | `c("Green","Red","ROMK")` | Output channel labels (input names are currently hardcoded as `"BFP"`, `"mCherry"`, `"GFP"`) |

**Known issue:** Bug #14 — input channel names (`"BFP"`, `"mCherry"`, `"GFP"`) are hardcoded; users with different DB channel names get all `NA`. A `channel_names` parameter needs to be added.

---

### `mergeSEandImg(se, df_img, tableType = "pa", selType, suffix = "hole")`

Joins cleaned imaging data into `colData` of an existing SCE.

| Argument | Default | Description |
|---|---|---|
| `se` | — | `SingleCellExperiment` |
| `df_img` | — | Output of `df_cleaned()` |
| `tableType` | `"pa"` | `"pa"` or `"coloc"` |
| `selType` | `c("Hole_ROI","background_ROI")` | Which ROI types to merge |
| `suffix` | `"hole"` | Suffix appended to channel names in colData (e.g. `"mCherry.hole"`) |

**Returns:** `se` with each imaging channel stored as a **nested DataFrame** in `colData` keyed as `"channelName.suffix"` (e.g. `colData(se)[["GFP.hole"]]`).

---

### `addThumbnailPaths(se, folder, ...)`

Scans a folder of JPG thumbnails, matches files to wells by `Plate_ID` + filename convention, and stores paths in `colData(se)$thumbnails` as a nested DataFrame.

**Filename convention expected:** `..._<Well>-<site>_<Channel>_<Class>_crop.jpg`
(e.g. `Exp_18T39265_A01-1_BF_class1_crop.jpg`)

**Known issue:** Bug #6 — older version overwrote `folder` with a hardcoded author path. Verify the current version uses the `folder` argument correctly.

---

### `plotThumbnails(se, wells, img_class)`

Renders a grid of BF/fluorescence thumbnail images for selected wells with contrast control.

---

### `mergeAnnotationCSV(se, csv, person = NULL)`

Merges an external annotation CSV into `colData`. Optionally filters to a specific annotator (`person` column).

---

## 4. Classification Pipeline

**File:** `R/classificationTools.R`

Replaces the deprecated `assign_cell_FINAL()` from `deprecated/Cellassignment_V2.R`.

### Pipeline overview

```
filterParticles()
    ↓
aggregateParticles()
    ↓
scoreParticles()
    ↓
classifyWells()
    ↓
mergeClassificationToSE()
```

Or use the all-in-one wrapper `classifyImgParticles()`.

---

### `filterParticles(df, method, threshold, group_vars, num_cols)`

Removes background / out-of-focus particles by intensity scaling within groups.

| Argument | Default | Description |
|---|---|---|
| `df` | — | Unaggregated particle data frame (output of `prepareImgDF`) |
| `method` | `"zscore"` | `"zscore"`: keep particles above group mean (threshold = 0); `"median_ratio"`: keep particles above `median * threshold` |
| `threshold` | `NULL` | Cut-off; `NULL` uses method default |
| `group_vars` | `c("Channel_Name","Plate_ID")` | Columns defining independent scaling groups |
| `num_cols` | `c("Area","Mean","IntDen","Number_of_Particles")` | Columns to set `NA` for rejected particles |

**Returns:** Input data frame with rejected particles' numeric columns set to `NA`; adds `Mean_scaled` column.

---

### `aggregateParticles(df, group_vars, agg_fun, scale_within, scale_center)`

Summarises filtered particles to one row per well × channel.

| Argument | Default | Description |
|---|---|---|
| `group_vars` | `c("Channel_Name","Plate_ID","Well")` | Aggregation grouping |
| `agg_fun` | `"mean"` | `"mean"`, `"median"`, or `"sum"` |
| `scale_within` | `"Channel_Name"` | Column(s) for optional within-group scaling of aggregated metrics |
| `scale_center` | `FALSE` | Whether to centre when scaling |

**Returns:** Data frame with columns `Mean_agg`, `Area_agg`, `normArea`, `n_particles` (and optionally `Mean_z`, `Area_z`, `normArea_z`). Wells where all particles were filtered score 0.

---

### `scoreParticles(agg_df, weights, score_group_vars, center)`

Rescales aggregated metrics within each channel group and combines them into a single `channel_score`.

| Argument | Default | Description |
|---|---|---|
| `weights` | `c(Mean=1, Area=1, normArea=1)` | Relative weights; normalised to sum to 1 internally |
| `score_group_vars` | `c("Channel_Name","Plate_ID")` | Independent scaling groups |
| `center` | `FALSE` | Use z-score (centred) vs uncentred (÷ sd) scaling. Centred is **not recommended** — it can score empty wells near-average |

**Returns:** `agg_df` with added `Mean_z`, `Area_z`, `normArea_z`, and `channel_score` columns.

---

### `classifyWells(score_df, delta, min_area, channel_labels)`

Pivots scores to wide format and calls each channel positive or negative per well.

| Argument | Default | Description |
|---|---|---|
| `delta` | `0.5` | Dominance tolerance: channel is positive if `score ≥ max_score - delta` |
| `min_area` | `0.1` | Minimum `normArea` (occupancy) for a positive call |
| `channel_labels` | `NULL` | Named vector mapping raw channel names to display labels (e.g. `c(C1="GFP", C2="mScarlett")`) |

**Returns:** Wide data frame (one row per well-plate) with `<channel>_score`, `<channel>_normArea`, `<channel>_positive`, `max_score`, and `Classification` columns (e.g. `"GFP+ mScarlett+"`, `"Negative"`, `"Multiple+"`).

---

### `classifyImgParticles(df, ...)` — All-in-one wrapper

Runs the full filter → aggregate → score → classify pipeline.

| Key arguments | Default | Description |
|---|---|---|
| `filter_method` | `"zscore"` | Passed to `filterParticles()` |
| `filter_threshold` | `NULL` | |
| `weights` | equal | Passed to `scoreParticles()` |
| `delta` | `0.5` | Passed to `classifyWells()` |
| `min_area` | `0.1` | |
| `channel_labels` | `NULL` | |
| `return_scores` | `FALSE` | Return scored long-format df instead of classification table |

---

### `mergeClassificationToSE(se, classification, col_name = NULL)`

Joins classification results to `colData` by `Well` + `Plate_ID`.

| Argument | Default | Description |
|---|---|---|
| `col_name` | `NULL` | `NULL` → add all new columns flat into colData; string → store as nested DataFrame under that name |

---

## 5. Clustering Pipeline

**File:** `R/clusteringTools.R`

Multi-modal unsupervised clustering combining electrophysiology and morphology
features. All functions work on `SingleCellExperiment` objects.

### Full pipeline

```
prepareClusterFeatures()
    ↓
optimalClusters()          ← optional; auto-selects k
    ↓
reducedDimMultimodal()     ← PCA + UMAP + t-SNE + clustering
    ↓
clusterHeatmap()           ← feature × well heatmap split by cluster
clusterSummary()           ← per-cluster boxplots + table
clusterMOFA()              ← optional MOFA2 multi-modal factors
```

Or use the all-in-one `clusterPipeline()`.

---

### `prepareClusterFeatures(se, ephys_cols, morpho_cols, assay_names, ephys_weight, morpho_weight, scale_method, na_action)`

Builds a scaled, modality-weighted feature matrix (wells × features).

| Argument | Default | Description |
|---|---|---|
| `ephys_cols` | `NULL` | colData scalar columns for ephys (e.g. `"Imax.minima"`, `"Vhalf.minima"`, `"Capacitance_mean"`) |
| `morpho_cols` | `NULL` | colData scalar columns for morphology (e.g. `"GFP_Mean_z"`, `"GFP_normArea"`) |
| `assay_names` | `NULL` | Assay names to include as full IV-curve features (one feature per sweep) |
| `ephys_weight` | `1` | Relative total-variance contribution of the ephys block |
| `morpho_weight` | `1` | Relative total-variance contribution of the morphology block |
| `scale_method` | `"zscore"` | `"zscore"`, `"minmax"`, or `"none"` — applied **within** each block |
| `na_action` | `"impute_median"` | `"impute_median"`, `"impute_zero"`, or `"omit"` |

**Returns:** Numeric matrix (wells × features). Ephys columns prefixed `e__`, morphology `m__`. Carries a `"modality"` attribute for downstream annotation.

**Modality weighting logic:** After scaling each block independently, the total variance of each block is rescaled to `weight × n_features`, so a block with `ephys_weight=2` contributes twice as much variance per feature to PCA as one with `morpho_weight=1`.

---

### `optimalClusters(feature_mat, k_range, methods, n_boot, seed)`

Runs k-means across a range of k and computes quality metrics.

| Argument | Default | Description |
|---|---|---|
| `k_range` | `2:8` | k values to evaluate |
| `methods` | `c("silhouette","wss")` | Metrics to compute; add `"gap"` explicitly (slow — runs `n_boot` bootstraps) |
| `n_boot` | `50` | Bootstrap replicates for gap statistic |

**Returns:** List with `$scores` (data.frame), `$plot` (ggplot with dashed lines at suggested k), `$suggested_k` (named integer vector, one per metric).

**Suggestion rules:**
- **silhouette** → k with highest average silhouette width
- **wss** → elbow via second derivative of within-cluster SS
- **gap** → Tibshirani firstSEmax rule

---

### `clusterSE(se, feature_mat, method, k, col_name, ...)`

Clusters wells and writes labels to `colData(se)[[col_name]]`.

| `method` | Requires `k` | Package needed |
|---|---|---|
| `"kmeans"` | Yes | base R |
| `"hierarchical"` | Yes | base R (`ward.D2` linkage) |
| `"louvain"` | No | `igraph` |
| `"gmm"` | Optional (NULL = BIC selection) | `mclust` (Suggests) |

---

### `reducedDimMultimodal(se, ..., method, k_clusters, cluster_method)`

Full replacement for the removed `reducedDim.Cellwise()`. Builds the modality-weighted feature matrix, runs PCA → t-SNE / UMAP, and calls `clusterSE()` on PCA coordinates.

**Adds to se:**
- `reducedDims(se)$PCA` — PCA scores matrix; carries `percentVar` attribute
- `reducedDims(se)$TSNE` — 2-column DataFrame
- `reducedDims(se)$UMAP` — 2-column DataFrame
- `colData(se)$cluster` — factor cluster labels

---

### `plotDimRed(se, redDim.method, colorColumns, point_size, alpha)`

Scatter plot of any `reducedDims` entry, coloured by one or more `colData` columns. Numeric columns use viridis (plasma), categorical use a discrete palette. Multiple `colorColumns` → `ggpubr::ggarrange` panel grid.

---

### `clusterHeatmap(se, feature_mat, cluster_col, scale, show_well_names, ...)`

`ComplexHeatmap` with features as rows and wells as columns, split by cluster.

- **Row annotation** — blue (`e__`) for ephys, orange (`m__`) for morphology (from `"modality"` attribute)
- **Column annotation** — cluster colour bar
- **`scale = TRUE`** — z-scores each feature row so colour reflects within-feature variation
- `e__` / `m__` prefixes are stripped for display

---

### `clusterSummary(se, feature_cols, cluster_col, plot_type, ncol)`

Per-cluster box / violin plots (one facet per feature) with jittered points, plus a summary table.

**Returns:** `list($plot, $table)` where `$table` has columns `cluster`, `feature`, `n`, `mean`, `sd`.

---

### `plotFeatureImportance(se, feature_mat, cluster_col, pcs, n_top)`

Two-panel feature importance analysis. Requires `reducedDimMultimodal()` or `clusterPipeline()` to have been run first (stores `PCA_model` in `metadata(se)`).

| Panel | What it shows | When it's high |
|---|---|---|
| **PCA loadings** | Signed contribution of each feature to PC1/PC2/PC3 | Feature drives global variance across all wells |
| **η² discriminability** | % of each feature's total variance that falls *between* clusters (SS_between / SS_total) | Feature cleanly separates the identified clusters |

A feature with **high loading + high η²** is the ideal cluster signature.
A feature with **high loading + low η²** varies a lot globally but doesn't separate your specific clusters.
A feature with **low loading + high η²** is quiet overall but perfectly sorts cells — often the most biologically specific.

```r
imp <- plotFeatureImportance(result$se, result$feature_mat,
                              pcs = 1:3, n_top = 15)
imp$combined_plot          # loadings (left) + discriminability (right)
imp$table                  # sorted by η², includes all PC loadings
```

**Returns:** `list($loadings_plot, $discriminability_plot, $combined_plot, $table)`.

---

### `clusterMOFA(se, ephys_cols, morpho_cols, assay_names, n_factors, scale_views, seed, ...)`

Optional MOFA2 wrapper. Builds separate ephys and morphology views (features × wells matrices), fits a multi-modal factor model, and stores factor scores in `reducedDims(se)[["MOFA"]]` and the model object in `metadata(se)[["MOFA_model"]]`.

Requires `BiocManager::install("MOFA2")` and a working Python environment. Factors loading on both modalities indicate shared cell-type signatures; modality-specific factors capture technical or biological variation unique to one modality.

---

### `clusterPipeline(se, ephys_cols, morpho_cols, ..., k, k_range, cluster_method, verbose)`

All-in-one wrapper. Runs steps 1–5 in sequence with progress messages.

```r
result <- clusterPipeline(
  se,
  ephys_cols    = c("Imax.minima", "Vhalf.minima", "Capacitance_mean"),
  morpho_cols   = c("GFP_Mean_z", "GFP_normArea",
                    "mCherry_Mean_z", "mCherry_normArea"),
  ephys_weight  = 1,
  morpho_weight = 1,
  k_range       = 2:6
)

result$optimal_k$plot       # silhouette + elbow curves
print(result$heatmap)       # feature heatmap
result$summary$plot         # per-cluster boxplots
result$summary$table        # per-cluster mean / SD
plotDimRed(result$se, "UMAP",
           colorColumns = c("cluster", "Imax.minima", "GFP_Mean_z"))
```

**Returns:** `list($se, $feature_mat, $optimal_k, $heatmap, $summary)`.

---

## 6. Utility & Analysis Functions

**File:** `R/utilityTools.R`

### `subsetSE(se, subset)`

Subsetting wrapper: evaluates a plain expression against `colData` columns (like `base::subset`). `NA`s are treated as `FALSE`.

```r
se_filtered <- subsetSE(se, Condition != "Empty" & QC == "Pass")
```

---

### `checkAssay(assayList, assayList.se)`

Validates that requested assay names exist in the SE; drops missing ones with a warning.

---

### `colAG(se, assayList, fun = mean, sweeps = row.names(se))`

Appends column-wise sweep aggregations to `colData`. For each assay in `assayList`, creates a `<assay>_mean` column in colData.

---

### `reducedDim.Cellwise(se, assayList, colNames, scaling, byRow, center, method, k_clusters)`

Computes PCA, t-SNE, and UMAP from selected assays and/or colData columns.

| Argument | Default | Description |
|---|---|---|
| `assayList` | `c()` | Assays to include as features |
| `colNames` | `c()` | colData columns to include as features |
| `scaling` | `"within"` | `"within"` scales per assay; `"global"` scales all features together |
| `center` | `FALSE` | Whether to centre when scaling |
| `method` | `c("pca","tsne","umap")` | Which reductions to compute |
| `k_clusters` | `3` | Number of k-means clusters (applied per reduction) |

**Returns:** `se` with `reducedDims` (PCA, TSNE, UMAP) and `cluster.umap`, `cluster.tsne`, `cluster.pca` in colData.

**Known issue:** TODO #10 — contains a `print(pca_data)` debug statement when `scaling == "global"`.

---

### `plotDimRed(se, redDim.method, colorColumns)`

Plots a dimensionality reduction result from `reducedDims(se)`.

| Argument | Default | Description |
|---|---|---|
| `redDim.method` | — | `"PCA"`, `"TSNE"`, or `"UMAP"` |
| `colorColumns` | `character()` | colData columns to use as additional colour panels |

**Returns:** `ggplot` or `ggpubr::ggarrange` of multiple coloured panels.

---

### `plotAssayVSSweeps(se, assayList, rowCol, colorGroup, wrapFormula, grouped)`

Plots assay values against sweep parameter (e.g. `V_Clamp`) as an IV-curve style figure.

| Argument | Default | Description |
|---|---|---|
| `assayList` | — | Assays to plot as y-axis |
| `rowCol` | — | rowData column to use as x-axis (e.g. `"V_Clamp"`) |
| `colorGroup` | `c()` | colData column to colour by |
| `grouped` | `TRUE` | `TRUE` = mean ± SEM per group; `FALSE` = individual well traces |
| `wrapFormula` | `NULL` | Formula for `facet_wrap()` |

---

### `get_metric(se, assay = "Minima", inward = TRUE)`

Extracts IV-curve metrics from a step-wise voltage clamp assay.

| Argument | Default | Description |
|---|---|---|
| `assay` | `"Minima"` | Name of assay to analyse |
| `inward` | `TRUE` | `TRUE` for inward currents (peak = min); `FALSE` for outward (peak = max) |

**Adds to colData:**
- `Imax.<assay>` — peak current amplitude
- `Vmax.<assay>` — holding potential at peak
- `Vhalf.<assay>` — half-activation voltage (spline-estimated)

**Uses** `sechm::meltSE` + `split()` + `lapply()` per Well+Plate_ID (bug #4 fixed).

---

### `fit_boltzmann_se(se, ...)`

Fits a Boltzmann sigmoidal curve to activation data. See roxygen docs for full parameter list.

---

## 7. Local Image Browser Module

**File:** `R/localImgBrowserModule.R`

Client-side Shiny module that lets users browse a local image folder from a **deployed** app (no server upload; images read entirely in the browser as blob: URLs).

### `localImgBrowserUI(id, label)`

Shiny UI function. Renders a button that opens the OS folder picker via `<input webkitdirectory>`.

### `localImgBrowserServer(id)`

Shiny server function. Returns a **reactive** with slots:
- `$n` — file count
- `$map` — named character vector: `"subdir|well.channel.class"` → `blob:URL`
- `$subdirs` — unique first-level subdirectory names (= Plate_IDs)
- `$channels` — unique channel names found
- `$classes` — unique class names found

### `filterImgMap(map, subdir, channel, img_class)`

Filters the `map` reactive by plate, channel, and/or image class. Returns a named vector suitable for `plateGridUI()`.

### `plateGridUI(url_map, coldata_colors, click_input_id)`

Renders an interactive plate grid (HTML) coloured by per-well metadata, with clickable wells.

**Filename convention:** `<well>.<channel>.<class>.jpg` (e.g. `A01.GFP.BF.jpg`), optionally under a subdirectory per plate.

---

## 8. Shiny App — tinySEV

**Files:** `R/ui.R`, `R/server.R`

### `tinySEV(objects, ...)`

Main entry point. Launches the interactive HT-APC dashboard.

```r
ephacRTools::tinySEV(objects = list("MyExperiment" = se))
```

### `tinySEV.ui(id)` / `tinySEV.server(id, objects, session)`

Module-level UI and server components for embedding `tinySEV` in a larger app.

**Dashboard sections:**

| Section | Features |
|---|---|
| **Prepare Object** | Overview, upload RDS/Excel, edit colData, sweep metadata |
| **Customize Object** | Define Conditions/Sweeps, rename/filter assays |
| **Plotting** | Plate heatmap, IV curve (sweep) plots |
| **Image Analysis** | Import SQLite DBs (`previewImgDB` → `prepareImgDF`), auto classification pipeline (6 steps), manual scoring, Plate Viewer with `localImgBrowserModule` |
| **Clustering** | PCA / t-SNE / UMAP via `reducedDim.Cellwise` |
| **Export** | Download SE (RDS), tables (CSV/Excel), plots |

**Optional features:** login (`shinyauthr` + `sodium`), memory monitoring, multi-dataset management.

---

## 9. Misc / Internal Helpers

**File:** `R/misc.R`

| Function | Description |
|---|---|
| `ag(df, cols, fun)` | Internal data aggregation helper. **Bug #3:** contains `print(c)` debug statement. |
| `dround(x, digits)` | Round to N significant digits. **Exported.** |
| `grepGene(pattern, names)` | Case-insensitive grep wrapper for gene/feature names. **Exported.** |
| `%\|\|%` | Null-coalescing operator (`x %\|\|% y` returns `y` if `x` is `NULL`). Internal. |

---

## 10. SCE Object Structure

```
SingleCellExperiment
├── assays            List of named matrices
│   ├── "Minima"      [sweeps × wells]  — per-sweep minimum current
│   ├── "Seal_Resistance"
│   ├── "Capacitance"
│   └── ... (one per Online Analysis output)
│
├── rowData           DataFrame [sweeps × metadata]
│   ├── Sweep         sweep number (integer)
│   ├── V_Clamp       holding potential (mV) — numeric or NA for ramp
│   └── Compound      condition label (collapsed or nested DataFrame)
│
├── colData           DataFrame [wells × metadata]
│   ├── Well          "A01" – "P24"
│   ├── Plate_ID      barcode from Nanion chip
│   ├── QC            QC flag from DataControl
│   ├── Row / Column  parsed from Well
│   │
│   ├── "GFP.hole"    nested DataFrame — imaging data per channel
│   │   ├── Mean, Area, IntDen, normArea
│   │   └── ... (from mergeSEandImg)
│   │
│   ├── thumbnails    nested DataFrame (from addThumbnailPaths)
│   │   └── <channel>.<class>  → file path
│   │
│   ├── Classification  "GFP+ mScarlett+", "Negative", …
│   └── Imax.minima, Vmax.minima, Vhalf.minima  (from get_metric)
│
├── reducedDims       List (after reducedDim.Cellwise)
│   ├── PCA           matrix [wells × components]
│   ├── TSNE          DataFrame [wells × 2]
│   └── UMAP          DataFrame [wells × 2]
│
└── colnames          interaction(Well, Plate_ID)  e.g. "A01.18T39265"
```

**Built-in datasets** (load with `data()`):
- `se_hAG` — human adrenal gland, voltage-ramp protocol, 3 imaging channels
- `se_iN` — iPSC-derived neurons, step-wise voltage clamp, 1 imaging channel
- `se_romk` — ROMK channel, used in the demo Shiny app

---

## 11. Known Issues

From `TODO.md` — open items as of 2026-04-01:

| # | Severity | File | Issue |
|---|---|---|---|
| #6 | 🔴 Bug | `imagingTools.R` | `addImgPaths()` may still reference hardcoded path — verify `folder` arg is used |
| #9 | 🔴 Bug | `prepareSE.R` | Non-unique colnames risk across plates (verify current fix is complete) |
| #3 | 🟡 Cleanup | `misc.R` | `print(c)` debug statement in `ag()` |
| #7 | 🟡 Cleanup | `imagingTools.R` | `addImgPaths()` → `imageval()` pipeline broken (returns NULL) |
| #8 | 🟡 Cleanup | `imagingTools.R` | Missing `EBImage::` / `grid::` namespace prefixes in `imageval()` |
| #10 | 🟡 Cleanup | multiple | Debug `print()`/`cat()` calls in `prepareDF`, `prepareSingleImgDF`, `reducedDim.Cellwise` |
| #11 | 🟡 Cleanup | `prepareSE.R` | `rm(list=setdiff(ls(),"new.df"))` inside function — dangerous, remove |
| #12 | 🔵 Infra | `DESCRIPTION` | `Description:` field empty, `License:` not set, `rlang` missing from Imports |
| #14 | 🟡 Cleanup | `imagingTools.R` | `df_cleaned()` hardcoded input channel names (`"BFP"`,`"mCherry"`,`"GFP"`) |
| #17 | 🔵 Infra | `vignettes/` | Vignette uses `devtools::load_all()` + relative paths — needs fixing for installed package |
| #18 | 🟡 Cleanup | `server.R` | `tryCatch` won't catch `shiny.silent.error` from `req()` |
| #19 | 🟡 Cleanup | `server.R` | Dead `.extract_plate_from_path()` function — remove |
