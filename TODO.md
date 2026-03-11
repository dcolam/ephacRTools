# ephacRTools — Bug & Improvement To-Do List

Generated from code review. Grouped by priority.

---

## 🔴 Fix Now — Correctness Bugs

### #1 · `group_assignment()` — 3 bugs
**File:** `R/RECOVERED_colDatTools.R`

- `pattern == "Conditions"` uses bare `df$Conditions` instead of `prepared_df$Conditions` → throws *"object 'df' not found"*
- `pattern == "Manual"` overwrites the whole data frame: `prepared_df <- manual_map[...]` → should be `prepared_df$Group <- manual_map[as.character(prepared_df$Well)]`
- Duplicate `else if (pattern == "Cycle")` block (~line 176) — dead code, remove it

---

### #2 · Duplicate `ag()` definition
**Files:** `R/misc.R` and `R/RECOVERED_colDatTools.R`

Identical function defined twice → R CMD CHECK error. Remove the copy in `RECOVERED_colDatTools.R`, keep the one in `misc.R`.

---

### #4 · `get_metric()` — plate loop only saves last plate's result
**File:** `R/utilityTools.R`

```r
for (i in seq_along(wells)) {
  for (p in seq_along(plate_ids)) {  # iterates all plates...
    ...
    results$Imax[i] <- Imax          # ...but saves only the last plate
  }
}
```

**Fix:** Export `get_metric_v2` as the canonical function (rename it `get_metric`), delete the broken original. `get_metric_v2` uses `split()` + `lapply()` correctly.

---

### #5 · `assign_cell_FINAL()` — misplaced parenthesis
**File:** `R/Cellassignment_V2.R`

```r
if (length(chan2analyze > 2))   # BUG: always TRUE
if (length(chan2analyze) > 2)   # correct
```

---

### #6 · `addImgPaths()` — ignores its parameter, uses hardcoded path
**File:** `R/imagingTools.R`

Function body immediately overwrites `parent_folder` with:
```r
pathImg <- "Y:\\ephacoffice\\DColameo\\DATA\\..."
```
Only works on author's machine. Replace `pathImg` with `parent_folder` and remove the hardcoded line.

---

### #9 · `prepareSE()` — non-unique column names across plates
**File:** `R/prepareSE.R`

```r
colnames(se) <- cd$Well   # "A01" appears in every plate → collision
```

`colData` rownames correctly use `interaction(Well, Plate_ID)` but are then overridden. Fix: keep colnames as `interaction(Well, Plate_ID)`.

---

## 🟡 Clean Up Soon

### #3 · `print(c)` debug statement in `ag()`
**File:** `R/misc.R` line ~9 — prints all numeric column names on every call. Remove it.

---

### #7 · Reconnect `addImgPaths()` → `imageval()` pipeline
**File:** `R/imagingTools.R`

`addImgPaths()` computes paths but never writes them back to the SE (returns `NULL` implicitly). `imageval()` then reads `se$Image_paths[idx]` which is always `NULL`.

**Fix:** Have `addImgPaths()` store results in `colData(se)$Image_paths` and return the updated SE.

---

### #8 · Missing package prefixes in `imageval()`
**File:** `R/imagingTools.R`

| Call | Should be |
|---|---|
| `readImage()` | `EBImage::readImage()` |
| `normalize()` | `EBImage::normalize()` |
| `rasterGrob()` | `grid::rasterGrob()` |
| `grid.arrange()` | `gridExtra::grid.arrange()` |

Also add `EBImage`, `grid`, `gridExtra` to `DESCRIPTION Imports`.

---

### #10 · Remove debug `cat()`/`print()` statements from production functions
Scattered across multiple files — inappropriate for a package:

| File | Statement |
|---|---|
| `R/prepareSE.R` `prepareDF()` | `cat("📦 Starting prepareDF")`, `cat("🧠 Memory...")`, `print(head(new.df, 3))` etc. |
| `R/prepareSE.R` `prepareMultipleDFs()` | `print("Excels Loaded")`, `print(safe_names)` |
| `R/utilityTools.R` `reducedDim.Cellwise()` | `print(pca_data)` when `scaling == "global"` |
| `R/utilityTools.R` `get_metric()` | `print(paste0("This Well failed:", well))` → replace with `warning()` |
| `R/imagingTools.R` `prepareSingleImgDF()` | `print(tbl)`, `print(names(tbl))`, `print("DB processed")` |
| `R/imagingTools.R` `df_cleaned()` | `print("DF cleaned")` |
| `R/imagingTools.R` `prepareImgDF()` | `print(df)` |

---

### #11 · Remove `rm(ls())` / `gc()` from `prepareDF()`
**File:** `R/prepareSE.R`

```r
rm(list = setdiff(ls(), "new.df"))
gc()
```

Calling `rm()` on `ls()` inside a function deletes all local variables mid-execution — dangerous and unnecessary. R garbage-collects locals on function exit automatically. Remove both lines.

---

### #13 · Remove `require()` calls from `plot_intensity_distribution()`
**File:** `R/Cellassignment_V2.R`

```r
require(ggplot2)   # wrong in package code — silently returns FALSE if missing
require(dplyr)
```

Both are already in `DESCRIPTION Imports`. Just remove these two lines.

---

### #14 · `df_cleaned()` — hardcoded input channel names
**File:** `R/imagingTools.R`

```r
ifelse(df$Channel_Name == "BFP",    channels[1],
  ifelse(df$Channel_Name == "mCherry", channels[2],
    ifelse(df$Channel_Name == "GFP",   channels[3], NA)))
```

The `channels` parameter only controls output labels; input names are hardcoded. Users with different database channel names get all `NA` silently.

**Fix:** Add a `channel_names = c("BFP", "mCherry", "GFP")` parameter that users can override.

---

### #15 · `filtered_df()` — double `Compound` assignment
**File:** `R/RECOVERED_colDatTools.R`

`Compound` is set unconditionally at line ~111, then immediately overwritten inside `if (ion == "Na")`. The first assignment is always a no-op. Remove the unconditional one.

---

### #16 · Delete empty placeholder files
**Files:** `R/cropTest.R`, `R/test.R`

Both files are empty. Remove from `R/`. If test scripts were intended, they belong in `tests/testthat/`.

---

### #19 · Remove dead `.extract_plate_from_path()` from `server.R`
**File:** `R/server.R` (~line 1906)

No longer called anywhere — replaced by `.ann_file_plate()`. Remove it.

---

## 🔵 Infrastructure

### #12 · Fix `DESCRIPTION`
**File:** `DESCRIPTION`

- [ ] Fill in the `Description:` field (currently *"More about what it does"*)
- [ ] Set a real `License:` (currently *"What license is it under?"*) — e.g. `MIT + file LICENSE`
- [ ] Remove `Remotes: dcolam/ephacRTools` (circular self-reference)
- [ ] Add missing `Imports`:
  - `tidyr` — `pivot_longer` / `pivot_wider` in `RECOVERED_colDatTools.R`
  - `jsonlite` — class buttons in `server.R`
  - `EBImage` — `imageval()`
  - `gridExtra` — `imageval()`
  - `scales` — `plot_well_assignment()`
  - `RColorBrewer` — `plot_well_assignment()`
  - `uwot` — `assign_cell_FINAL()` clustering mode
  - `rlang` — `%||%` used throughout `server.R`

---

### #17 · Fix vignette so it works after package install
**File:** `vignettes/getting-started.Rmd`

- Replace `devtools::load_all()` with `library(ephacRTools)`
- Replace relative paths like `"../data-raw/iNeurons/"` with bundled sample datasets (`data(se_iN)`)

---

### #18 · Fix `ann_plate_dirs()` — `tryCatch` won't catch `req()` failures
**File:** `R/server.R`

`req()` throws a `shiny.silent.error`, not a base error, so `tryCatch(..., error=function(e){})` never triggers the fallback when no SE is loaded.

**Fix:**
```r
tryCatch(img_plate_dirs(),
  error              = function(e) list(),
  shiny.silent.error = function(e) list())
```
Or check `isTruthy(SE())` before calling `img_plate_dirs()`.

---

## 📋 Summary Table

| # | File | Severity | Topic |
|---|---|---|---|
| 1 | `RECOVERED_colDatTools.R` | 🔴 Bug | `group_assignment()` — 3 bugs |
| 2 | `misc.R` / `RECOVERED_colDatTools.R` | 🔴 Bug | Duplicate `ag()` |
| 3 | `misc.R` | 🟡 Cleanup | `print(c)` debug statement |
| 4 | `utilityTools.R` | 🔴 Bug | `get_metric()` plate loop |
| 5 | `Cellassignment_V2.R` | 🔴 Bug | `length(x > 2)` parenthesis |
| 6 | `imagingTools.R` | 🔴 Bug | `addImgPaths()` hardcoded path |
| 7 | `imagingTools.R` | 🟡 Cleanup | `addImgPaths()` → `imageval()` broken pipeline |
| 8 | `imagingTools.R` | 🟡 Cleanup | Missing `EBImage::`/`grid::` prefixes |
| 9 | `prepareSE.R` | 🔴 Bug | Non-unique colnames across plates |
| 10 | multiple | 🟡 Cleanup | Remove debug `print()`/`cat()` calls |
| 11 | `prepareSE.R` | 🟡 Cleanup | `rm(ls())` inside function |
| 12 | `DESCRIPTION` | 🔵 Infra | License, Imports, circular Remotes |
| 13 | `Cellassignment_V2.R` | 🟡 Cleanup | `require()` in package code |
| 14 | `imagingTools.R` | 🟡 Cleanup | `df_cleaned()` hardcoded channel names |
| 15 | `RECOVERED_colDatTools.R` | 🟡 Cleanup | `filtered_df()` double assignment |
| 16 | `R/cropTest.R`, `R/test.R` | 🟡 Cleanup | Empty placeholder files |
| 17 | `vignettes/` | 🔵 Infra | Vignette not installable |
| 18 | `server.R` | 🟡 Cleanup | `tryCatch` won't catch `req()` |
| 19 | `server.R` | 🟡 Cleanup | Dead `.extract_plate_from_path()` |
