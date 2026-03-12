# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

`ephacRTools` v0.3.0 is an R package for analyzing **High-Throughput Automated Patch Clamp (HT-APC)** data from Syncropatch/Nanion instruments. It:

1. Converts DataControl Excel exports (sweep-wise electrophysiology data) into `SingleCellExperiment` objects
2. Integrates imaging data from Cluster Analysis SQLite databases
3. Provides a particle-level → well-level classification pipeline for cell typing
4. Offers interactive visualization via a `tinySEV` Shiny dashboard
5. Supports deployment to Posit Connect Cloud

**Authors**: David Colameo (aut, cre), Athena Schumacher (aut), Jiayi Wang (aut)

## Development Commands

```r
# Load package during development (preferred over library())
devtools::load_all()

# Build documentation from roxygen2 comments
devtools::document()

# Run checks
devtools::check()

# Run tests
devtools::test()

# Build and install
devtools::install()

# Build vignette
devtools::build_vignettes()

# Launch the Shiny app
source("app/app.R")
# Or directly:
ephacRTools::tinySEV(objects = list("MyData" = se_object))
```

## Architecture

### Data Flow

```
Syncropatch DataControl Excel (.xlsx)
    → prepareDF() / prepareMultipleDFs()   [R/prepareSE.R]
    → prepareSE()                           [R/prepareSE.R]
    → SingleCellExperiment (SCE) object

SQLite DB (Cluster Analysis output)
    → prepareSingleImgDF() / prepareImgDF() [R/imagingTools.R]
    → df_cleaned()                          [R/imagingTools.R]
    → mergeSEandImg()                       [R/imagingTools.R]
    → SCE with nested imaging DataFrames in colData

JPG thumbnails folder
    → addThumbnailPaths()                   [R/imagingTools.R]
    → SCE colData$thumbnails (nested DataFrame with file paths)
    → plotThumbnails()                      [R/imagingTools.R]  ← grid display

Particle imaging data (from prepareImgDF)
    → filterParticles()                     [R/classificationTools.R]
    → aggregateParticles()                  [R/classificationTools.R]
    → scoreParticles()                      [R/classificationTools.R]
    → classifyWells()                       [R/classificationTools.R]
    → mergeClassificationToSE()             [R/classificationTools.R]
    (or all-in-one: classifyImgParticles()) [R/classificationTools.R]

SCE object
    → colAG(), reducedDim.Cellwise(), get_metric()  [R/utilityTools.R]
    → fit_boltzmann_se()                            [R/utilityTools.R]
    → tinySEV() Shiny app                           [R/ui.R + R/server.R]
```

### SCE Object Structure

- **assays**: Named matrices; each assay = one electrophysiology measurement (e.g., "Minima", "Seal_Resistance", Capacitance), rows = sweeps, columns = wells
- **rowData**: Sweep-level metadata (`sweep`, `V_Clamp` voltage, compound info). When compound layout varies across wells, stored as nested DataFrames.
- **colData**: Well-level metadata (`Well`, `Plate_ID`, `QC`, `Row`, `Column`). Imaging data stored as **nested DataFrames** per channel (e.g., `colData(se)[["mCherry.hole"]]` = DataFrame with `Area`, `Mean`, `IntDen`, `normArea` columns). Classification results added as flat columns or nested DataFrame.

### Key Files

| File | Purpose |
|------|---------|
| `R/prepareSE.R` | Parse DataControl Excel → long-format df → SCE. Entry: `prepareSE()` |
| `R/imagingTools.R` | Parse SQLite imaging DBs, merge into SCE colData; thumbnail path tracking. Entry: `prepareImgDF()` + `mergeSEandImg()` + `addThumbnailPaths()` |
| `R/classificationTools.R` | **NEW** Modular cell-type classification pipeline: filter → aggregate → score → classify → merge. Entry: `classifyImgParticles()` |
| `R/localImgBrowserModule.R` | **NEW** Client-side Shiny module for browser-based folder selection (works on deployed servers). Exports: `localImgBrowserUI()`, `localImgBrowserServer()`, `filterImgMap()`, `plateGridUI()` |
| `R/utilityTools.R` | Aggregation (`colAG`), dim reduction (`reducedDim.Cellwise`), IV-curve metrics (`get_metric`, `fit_boltzmann_se`), plotting (`plotDimRed`, `plotAssayVSSweeps`) |
| `R/ui.R` / `R/server.R` | `tinySEV` Shiny dashboard — interactive viewer with heatmaps, IV plots, image browser, classification pipeline, clustering, export |
| `R/misc.R` | Internal helpers: `ag()`, `dround()`, `grepGene()`, `%||%` |
| `R/data.R` | Dataset documentation for `se_hAG`, `se_iN` |
| `R/help.R` | Internal: modal help dialogs for the Shiny app |
| `app/app.R` | Standalone Shiny app launcher using built-in `se_romk` dataset |
| `app/imgBrowserDemo.R` | **NEW** Standalone demo app for `localImgBrowserUI` module |
| `deploy.R` | **NEW** Deployment script for Posit Connect Cloud (renv, manifest, rsconnect) |
| `deprecated/` | Old code (e.g., `Cellassignment_V2.R` with `assign_cell_FINAL()`) — superseded by `classificationTools.R` |

### Sample Datasets

Three built-in datasets (loaded with `data()`):
- `se_hAG`: Human adrenal gland, voltage-ramp protocol, 4 measurements + 3 imaging datasets
- `se_iN`: iPSC-derived neurons, step-wise voltage clamp (-80mV to +20mV), 1 imaging dataset
- `se_romk`: Used in the demo Shiny app (`app/app.R`)

## Classification Pipeline Details (`R/classificationTools.R`)

New modular pipeline replacing the older `assign_cell_FINAL()` approach:

```r
# All-in-one
result <- classifyImgParticles(df,
  filter_method = "zscore",    # or "median_ratio"
  filter_threshold = 0,         # zscore: keep > mean; median_ratio: keep > 1/3 of median
  weights = c(Mean=1, Area=1, normArea=1),
  delta = 0.5,                  # score within delta of max → positive
  min_area = 0.1,               # minimum normArea to be positive
  channel_labels = c(C1="GFP", C2="mCherry")
)
se <- mergeClassificationToSE(se, result)

# Step-by-step
filtered <- filterParticles(df)        # remove background particles
agg      <- aggregateParticles(filtered) # per well/channel summary
scored   <- scoreParticles(agg)         # uncentred scaling + weighted composite
classes  <- classifyWells(scored)       # pivot wide, threshold → "GFP+ mCherry+"
se       <- mergeClassificationToSE(se, classes)
```

**Scoring logic**: Each metric (Mean, Area, normArea) is scaled 0–1 within channel/plate group → combined as weighted sum (`channel_score`). A channel is "positive" if `score ≥ (max_score - delta)` AND `normArea > min_area`.

## Imaging Pipeline Details (`R/imagingTools.R`)

- Imaging data from **Cluster Analysis** software, exported as SQLite databases
- `prepareSingleImgDF()`: queries `Particle_Analysis_Table` (analysis="pa") or `Coloc_Analysis_Table` (analysis="coloc") joined with measurement tables; computes Mean, Area, IntDen, normArea per well/channel
- `df_cleaned()`: labels ROIs as `"Hole_ROI"` (smallest area) vs `"background_ROI"` (largest area) per selection
- After `mergeSEandImg()`: each imaging channel becomes a nested DataFrame in colData keyed as `"channelName.suffix"` (e.g., `"mCherry.hole"`)
- `addThumbnailPaths(se, folder)`: scans for JPGs, matches to wells by Plate_ID + filename convention, stores paths in `colData(se)$thumbnails` as nested DataFrame
- `plotThumbnails(se, wells, img_class)`: renders BF/fluorescence grids with contrast control

## Local Image Browser Module (`R/localImgBrowserModule.R`)

Client-side file browser for Shiny apps — reads files **entirely in the browser** via JavaScript (no server upload). Works on deployed Posit Connect apps.

```r
# In UI:
localImgBrowserUI("imgBrowser")

# In server:
img <- localImgBrowserServer("imgBrowser")
# Returns reactive: img$n (count), img$map (named vec "subdir|well.channel.class" → blob: URL),
#                   img$subdirs, img$channels, img$classes

# Filter and render plate grid:
url_map <- filterImgMap(img$map(), subdir="Plate1", channel="GFP", img_class="BF")
plateGridUI(url_map, coldata_colors = color_vec, click_input_id = "clicked_well")
```

- Uses `<input type="file" webkitdirectory>` triggered via JS button
- First-level subdirectory = Plate_ID
- Filename convention: `well.channel.class.jpg` (e.g., `A01.GFP.BF.jpg`)

## tinySEV Shiny App

Dashboard with the following sections (sidebar):
1. **Prepare Object** → Overview, Upload RDS/Excel, Column Data editor, Sweep metadata
2. **Customize Object** → Define Conditions, Define Sweeps, Change Assays, Filtering
3. **Plotting** → Plate Heatmap, IV Curve (Sweep) Plots
4. **Image Analysis** → Import Data (SQLite/JPG), Auto Classification (6-step pipeline), Manual Scoring, Plate Viewer
5. **Clustering** → PCA/t-SNE/UMAP dimensionality reduction
6. **Export** → Download SE, tables, plots

Features: optional login (shinyauthr + sodium), memory monitoring, multi-SE dataset management, client-side image browsing.

## Testing

```r
tests/
├── testthat.R
└── testthat/
    ├── test-filterImgMap.R   # 10 tests for filterImgMap (channel/class/subdir filtering, dedup)
    ├── test-misc.R           # dround, grepGene, etc.
    └── test-prepareSE.R      # prepareSE pipeline
```

## Deployment

See `deploy.R` for Posit Connect Cloud workflow:
1. Rebuild docs + reinstall from GitHub
2. Snapshot `renv.lock` from `app/` directory
3. Write `manifest.json` (then remove BiocVersion entry)
4. Deploy via `rsconnect::deployApp()`

## Document Generation

All exported functions use roxygen2 (`#'`) comments. Run `devtools::document()` after any signature or documentation change. The `NAMESPACE` file is auto-generated — never edit manually.

## Exported Functions Reference

| Category | Function | File |
|----------|----------|------|
| Data prep | `prepareDF`, `prepareMultipleDFs`, `prepareSE` | prepareSE.R |
| Imaging | `prepareSingleImgDF`, `prepareImgDF`, `df_cleaned`, `mergeSEandImg` | imagingTools.R |
| Imaging | `addThumbnailPaths`, `plotThumbnails` | imagingTools.R |
| Classification | `filterParticles`, `aggregateParticles`, `scoreParticles` | classificationTools.R |
| Classification | `classifyWells`, `classifyImgParticles`, `mergeClassificationToSE` | classificationTools.R |
| Image browser | `localImgBrowserUI`, `localImgBrowserServer`, `filterImgMap`, `plateGridUI` | localImgBrowserModule.R |
| Utilities | `checkAssay`, `colAG`, `reducedDim.Cellwise`, `get_metric`, `fit_boltzmann_se` | utilityTools.R |
| Utilities | `plotDimRed`, `plotAssayVSSweeps`, `dround`, `grepGene` | utilityTools.R / misc.R |
| Shiny | `tinySEV`, `tinySEV.ui`, `tinySEV.server` | ui.R / server.R |
