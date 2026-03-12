# Image Classification Pipeline

`R/classificationTools.R` — particle-level cell classification from Cluster Analysis imaging data.

## Overview

The pipeline classifies each well as positive or negative for one or more fluorescent channels using particle-level data exported from the Cluster Analysis SQLite databases. It works in four steps:

```
prepareImgDF()  ──►  filterParticles()  ──►  aggregateParticles()
                                                      │
                                                      ▼
mergeClassificationToSE()  ◄──  classifyWells()  ◄──  scoreParticles()
```

The high-level wrapper `classifyImgParticles()` runs all four middle steps in one call.

---

## Functions

### `filterParticles()`

**Purpose:** Remove background / out-of-focus particles before aggregation by blanking their numeric columns.

**Signature:**
```r
filterParticles(
  df,
  method     = c("zscore", "median_ratio"),
  threshold  = NULL,
  group_vars = c("Channel_Name", "Plate_ID"),
  num_cols   = c("Area", "Mean", "IntDen", "Number_of_Particles")
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `df` | data.frame | — | Unaggregated particle table. Must contain `Mean` and the columns in `group_vars`. |
| `method` | character | `"zscore"` | Scaling strategy. See details below. |
| `threshold` | numeric or NULL | `NULL` | Rejection cut-off. `NULL` uses the method's default. |
| `group_vars` | character vector | `c("Channel_Name", "Plate_ID")` | Columns used to define groups for scaling. |
| `num_cols` | character vector | `c("Area", "Mean", "IntDen", "Number_of_Particles")` | Columns set to `NA` for rejected particles. |

**Methods:**
- **`"zscore"`** (default): standardises `Mean` (mean=0, sd=1) within each group. Particles with `Mean_scaled < threshold` are rejected. Default threshold = `0` (below-average particles are dropped).
- **`"median_ratio"`**: computes `Mean / median(Mean)` within each group. Particles below `threshold` are rejected. Default threshold = `1/3` (particles at less than one-third of the group median are dropped).

**Returns:** The input data frame with `num_cols` set to `NA` for rejected rows, plus an appended `Mean_scaled` column.

**Example:**
```r
df_filtered <- filterParticles(df_raw, method = "zscore", threshold = 0)

# Less aggressive: keep particles above 1/4 of the median
df_filtered <- filterParticles(df_raw, method = "median_ratio", threshold = 0.25)
```

---

### `aggregateParticles()`

**Purpose:** Collapse the particle-level table to one row per well and channel, computing summary statistics used for scoring.

**Signature:**
```r
aggregateParticles(
  df,
  group_vars = c("Channel_Name", "Plate_ID", "Well")
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `df` | data.frame | — | Particle table, typically the output of `filterParticles()`. |
| `group_vars` | character vector | `c("Channel_Name", "Plate_ID", "Well")` | Columns defining aggregation groups. |

**Behaviour:**
- If a `CorrSel` column is present, only rows with `CorrSel == "Hole_ROI"` are used (the well ROI, not the surrounding background).
- All well × channel groups present in the input are preserved in the output. Groups where every particle was rejected by `filterParticles()` (all `Mean == NA`) appear as rows with `NA` stats. This ensures they score 0 in `scoreParticles()` and are ultimately classified as `"Negative"` rather than silently dropped.

**Output columns:**

| Column | Description |
|--------|-------------|
| `Mean_agg` | Mean fluorescence intensity of accepted particles. |
| `Area_agg` | Mean particle area (pixels²). |
| `normArea` | Occupancy: `sum(Area) / mean(Selection_Area)`. Values near 1 mean particles cover the entire hole ROI. |
| `n_particles` | Count of accepted (non-NA) particles. |

**Example:**
```r
agg <- aggregateParticles(df_filtered)
# or with custom grouping:
agg <- aggregateParticles(df_filtered, group_vars = c("Channel_Name", "Well"))
```

---

### `scoreParticles()`

**Purpose:** Convert aggregated per-channel metrics into a single composite `channel_score` per well, enabling cross-channel comparison.

**Signature:**
```r
scoreParticles(
  agg_df,
  weights = c(Mean = 1, Area = 1, normArea = 1)
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `agg_df` | data.frame | — | Output of `aggregateParticles()`. |
| `weights` | named numeric | `c(Mean=1, Area=1, normArea=1)` | Relative importance of each metric. Normalised internally to sum to 1. |

**How it works:**
1. Each metric (`Mean_agg`, `Area_agg`, `normArea`) is scaled **without centering** (divided by its SD) **within each channel group** (default: per `Channel_Name × Plate_ID`). Scaling per-channel means a dim channel and a bright channel each have scores relative to their own well populations, making cross-channel comparisons inside `classifyWells()` fair.
2. `NA` values (wells where all particles were filtered) are replaced with `0` before combining.
3. The weighted sum is computed: `channel_score = w_Mean * Mean_z + w_Area * Area_z + w_normArea * normArea_z`.

**Output columns added:** `Mean_z`, `Area_z`, `normArea_z`, `channel_score`.

**Tuning advice:**
- If intensity is the most reliable signal, upweight `Mean`: `c(Mean = 2, Area = 1, normArea = 1)`.
- If you only care about coverage (occupancy), use: `c(Mean = 0, Area = 0, normArea = 1)`.

**Example:**
```r
# Equal weights (default)
scored <- scoreParticles(agg)

# Favour intensity over area
scored <- scoreParticles(agg, weights = c(Mean = 3, Area = 1, normArea = 2))
```

---

### `classifyWells()`

**Purpose:** Turn per-channel scores into binary positive/negative calls and a final combined `Classification` label per well.

**Signature:**
```r
classifyWells(
  score_df,
  delta          = 0.5,
  min_area       = 0.1,
  channel_labels = NULL
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `score_df` | data.frame | — | Output of `scoreParticles()`. |
| `delta` | numeric | `0.5` | Dominance tolerance. A channel is called positive only if its score is within `delta` of the well's maximum score across all channels. |
| `min_area` | numeric | `0.1` | Minimum occupancy (`normArea`) required for a positive call. Prevents low-area noise from being classified as positive. |
| `channel_labels` | named character or NULL | `NULL` | Rename channels in the output. E.g. `c(C1 = "GFP", C2 = "mScarlett")`. |

**How positivity is determined:**

For each well, the maximum channel score is found (`max_score = max(all channel scores)`). Channel `ch` is positive if:

```
score_ch  >=  max_score - delta   AND   normArea_ch  >  min_area
```

- `delta` controls how similar a secondary channel's score needs to be to the dominant channel before it is also called positive. A large `delta` (e.g. `1.0`) makes calling co-positives more permissive; a small `delta` (e.g. `0.1`) only calls the single best-scoring channel positive.
- `min_area` is an absolute floor that prevents faint/sparse signal from being called positive regardless of score.

**Output columns:**

| Column | Description |
|--------|-------------|
| `Well`, `Plate_ID` | Identifiers. |
| `<channel>_score` | Composite score per channel. |
| `<channel>_normArea` | Occupancy per channel. |
| `<channel>_positive` | `TRUE`/`FALSE` per channel. |
| `max_score` | Highest score in that well. |
| `Classification` | Label string, e.g. `"GFP+"`, `"GFP+ mScarlett+"`, `"Negative"`, `"Multiple+"`. |

**Classification label rules:**
- No channel positive → `"Negative"`
- All channels positive (>2 channels) → `"Multiple+"`
- Otherwise → channel names of positive channels joined with `" "`, each suffixed with `"+"`, e.g. `"GFP+ mScarlett+"`

**Example:**
```r
cls <- classifyWells(scored, delta = 0.5, min_area = 0.1,
                     channel_labels = c(C1 = "GFP", C2 = "mScarlett"))
table(cls$Classification)
```

---

### `classifyImgParticles()` _(main wrapper)_

**Purpose:** Convenience function that chains `filterParticles → aggregateParticles → scoreParticles → classifyWells` in a single call.

**Signature:**
```r
classifyImgParticles(
  df,
  filter_method     = c("zscore", "median_ratio"),
  filter_threshold  = NULL,
  filter_group_vars = c("Channel_Name", "Plate_ID"),
  agg_group_vars    = c("Channel_Name", "Plate_ID", "Well"),
  weights           = c(Mean = 1, Area = 1, normArea = 1),
  delta             = 0.5,
  min_area          = 0.1,
  channel_labels    = NULL,
  return_scores     = FALSE
)
```

**Parameters:**

| Parameter | Passed to | Description |
|-----------|-----------|-------------|
| `df` | `filterParticles` | Unaggregated particle data frame. Must be filtered to `Image_Type == "fluor"` beforehand. |
| `filter_method` | `filterParticles` | `"zscore"` or `"median_ratio"`. |
| `filter_threshold` | `filterParticles` | Cut-off value; `NULL` uses method default. |
| `filter_group_vars` | `filterParticles` | Groups for scaling. |
| `agg_group_vars` | `aggregateParticles` | Groups for aggregation. |
| `weights` | `scoreParticles` | Metric weights. |
| `delta` | `classifyWells` | Dominance tolerance. |
| `min_area` | `classifyWells` | Occupancy floor. |
| `channel_labels` | `classifyWells` | Channel rename map. |
| `return_scores` | — | If `TRUE`, return the scored long-format table instead of the final classification. Useful for inspection/debugging. |

**Returns:** Wide-format classification data frame (default) or scored long-format data frame (`return_scores = TRUE`).

**Typical usage:**
```r
# Step 1: Load raw particle data (unaggregated, fluor only)
df_raw <- prepareImgDF(
  pathDB     = "path/to/db.sqlite",
  aggregate  = FALSE,
  cleanNames = FALSE
)
df_raw <- subset(df_raw, Image_Type == "fluor")

# Step 2: Classify
cls <- classifyImgParticles(
  df_raw,
  filter_method  = "zscore",
  delta          = 0.5,
  min_area       = 0.1,
  channel_labels = c(C1 = "GFP", C2 = "mScarlett")
)

# Inspect
table(cls$Classification)
#> GFP+   mScarlett+   GFP+ mScarlett+   Negative
#>   42           18               31         153

# Step 3: Merge into SE
se <- mergeClassificationToSE(se, cls)
```

**Debugging intermediate steps:**
```r
# Inspect scores before classification
scored <- classifyImgParticles(df_raw, return_scores = TRUE)
hist(scored$channel_score)

# Run step by step
df_filt   <- filterParticles(df_raw)
df_agg    <- aggregateParticles(df_filt)
df_scored <- scoreParticles(df_agg)
cls       <- classifyWells(df_scored, delta = 0.5, min_area = 0.1)
```

---

### `mergeClassificationToSE()`

**Purpose:** Join classification results into the `colData` of a `SummarizedExperiment` object.

**Signature:**
```r
mergeClassificationToSE(
  se,
  classification,
  col_name = NULL
)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `se` | SummarizedExperiment | — | SE object with `Well` and `Plate_ID` in `colData`. |
| `classification` | data.frame | — | Output of `classifyImgParticles()` or `classifyWells()`. |
| `col_name` | character or NULL | `NULL` | If `NULL`, all new columns are added flat to `colData`. If a string, all new columns are stored as a nested `DataFrame` under that name. |

**Flat merge (default):**
```r
se <- mergeClassificationToSE(se, cls)
# Access directly:
colData(se)$Classification
colData(se)$GFP_score
colData(se)$GFP_positive
```

**Nested merge:**
```r
se <- mergeClassificationToSE(se, cls, col_name = "img_classification")
# Access via nested DataFrame:
colData(se)[["img_classification"]]$Classification
```

---

## Parameter Tuning Guide

### Choosing a filter method

| Scenario | Recommended method | Notes |
|----------|--------------------|-------|
| Most wells have cells | `"zscore"`, threshold `0` | Drops below-average particles; robust to bright outliers |
| Sparse staining, many empty wells | `"median_ratio"`, threshold `0.25–0.5` | Less aggressive; median is 0 in empty wells so effectively skips them |
| Very noisy background | `"zscore"`, threshold `0.5–1.0` | Stricter; only keeps clearly above-average particles |

### Choosing `delta`

| `delta` | Effect |
|---------|--------|
| `0.1` | Only the single highest-scoring channel per well is called positive (strict single-positive) |
| `0.5` (default) | Channels within half a SD-unit of the top score can also be positive (tolerates similar-scoring doublets) |
| `1.0+` | Very permissive; most channels in a well will be co-positive |

### Choosing `min_area`

`min_area` is in the same units as `normArea` (fraction of selection area). A value of `0.1` means at least 10% of the hole ROI must be covered by accepted particles for a channel to be called positive. Increase this if you are seeing false positives from a few sparse bright particles; decrease it if true positive cells are small.
