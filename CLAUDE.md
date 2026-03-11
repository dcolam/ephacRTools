# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Package Does

`ephacRTools` is an R package for analyzing **High-Throughput Automated Patch Clamp (HT-APC)** data from Syncropatch/Nanion instruments. It converts DataControl Excel exports (sweep-wise electrophysiology data) into `SingleCellExperiment` objects (built on `SummarizedExperiment`) and provides tools for quality control, cell assignment via imaging, dimensionality reduction, and interactive visualization.

## Development Commands

```r
# Load package during development (preferred over library())
devtools::load_all()

# Build documentation from roxygen2 comments
devtools::document()

# Run checks
devtools::check()

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

SCE object
    → assign_cell_FINAL() / assign_cell_FINAL_valcol()  [R/Cellassignment_V2.R]
    → colAG(), reducedDim.Cellwise(), get_metric()       [R/utilityTools.R]
    → tinySEV() Shiny app                               [R/ui.R + R/server.R]
```

### SCE Object Structure

- **assays**: Named matrices where each assay = one electrophysiology measurement (e.g., "Minima", "Seal_Resistance", capacitance, etc.), rows = sweeps, columns = wells
- **rowData**: Sweep-level metadata (sweep number, `V_Clamp` voltage, compound info). When compound layout varies across wells, stored as nested DataFrames.
- **colData**: Well-level metadata (`Well`, `Plate_ID`, `QC`, `Row`, `Column`). Imaging data is stored as **nested DataFrames** per channel (e.g., `colData(se)[["mCherry"]]` is a DataFrame with `Area`, `Mean`, `IntDen` columns).

### Key Files

| File | Purpose |
|------|---------|
| `R/prepareSE.R` | Parse DataControl Excel → long-format df → SCE. Entry point: `prepareSE()` |
| `R/imagingTools.R` | Parse SQLite imaging DBs, merge into SCE colData. Entry point: `prepareImgDF()` + `mergeSEandImg()` |
| `R/Cellassignment_V2.R` | Assign cell types (green/red/black) from imaging channels via manual thresholds or k-means/UMAP clustering |
| `R/utilityTools.R` | Aggregation (`colAG`), dimensionality reduction (`reducedDim.Cellwise`), IV-curve metrics (`get_metric`, `get_metric_v2`, `fit_boltzmann_se`), plotting (`plotDimRed`, `plotAssayVSSweeps`) |
| `R/ui.R` / `R/server.R` | `tinySEV` Shiny app — interactive viewer for SCE objects with heatmaps, plots, and file download |
| `R/misc.R` | Internal helpers: `ag()` (aggregate), `dround()`, `grepGene()` |
| `R/RECOVERED_colDatTools.R` | Older/recovered utilities: `prep_excel()`, `prep_df()` — legacy approach superseded by `prepareSE()` |
| `R/data.R` | Dataset documentation for built-in samples `se_hAG` and `se_iN` |
| `app/app.R` | Standalone Shiny app launcher using built-in `se_romk` dataset |

### Sample Datasets

Two built-in datasets (loaded with `data()`):
- `se_hAG`: Human adrenal gland, voltage-ramp protocol, 4 measurements + 3 imaging datasets
- `se_iN`: iPSC-derived neurons, step-wise voltage clamp (-80mV to +20mV), 1 imaging dataset
- `se_romk`: Used in the demo Shiny app

### Imaging Pipeline Details

- Imaging data comes from **Cluster Analysis** software, exported as SQLite databases
- `prepareSingleImgDF()` queries either `Particle_Analysis_Table` (particle analysis) or `Coloc_Analysis_Table` (colocalization) joined with measurement tables
- `df_cleaned()` labels ROIs as `"Hole_ROI"` (smallest area) vs `"background_ROI"` (largest area) per selection
- After `mergeSEandImg()`, each imaging channel becomes a nested DataFrame in colData, keyed by channel name (e.g., `"mCherry.hole"`)
- Cell assignment in `assign_cell_FINAL()` uses log1p-scaled per-channel Mean intensities with either manual thresholds (producing "green"/"red"/"black" labels) or UMAP + k-means clustering

### Document Generation

All exported functions use roxygen2 (`#'`) comments. Run `devtools::document()` after any signature or documentation change. The `NAMESPACE` file is auto-generated — do not edit manually.
