library(testthat)
library(ephacRTools)
library(SummarizedExperiment)
library(SingleCellExperiment)

# ── SCE object structure ─────────────────────────────────────────────────────

test_that("se_iN is a SingleCellExperiment", {
  data(se_iN)
  expect_s4_class(se_iN, "SingleCellExperiment")
})

test_that("se_hAG is a SingleCellExperiment", {
  data(se_hAG)
  expect_s4_class(se_hAG, "SingleCellExperiment")
})

test_that("required colData columns are present", {
  data(se_iN)
  cd <- colnames(colData(se_iN))
  expect_true("Well"     %in% cd)
  expect_true("Plate_ID" %in% cd)
  expect_true("QC"       %in% cd)
})

test_that("column names are unique well.plate_id composites", {
  data(se_iN)
  ids <- colnames(se_iN)
  expect_equal(length(ids), length(unique(ids)))
  # Composites contain a "." separator
  expect_true(all(grepl("\\.", ids)))
})

test_that("column names match interaction(Well, Plate_ID)", {
  data(se_iN)
  cd       <- as.data.frame(colData(se_iN))
  expected <- as.character(interaction(cd$Well, cd$Plate_ID))
  expect_equal(colnames(se_iN), expected)
})

test_that("assay matrices have consistent dimensions", {
  data(se_iN)
  for (a in assayNames(se_iN)) {
    m <- assay(se_iN, a)
    expect_equal(ncol(m), ncol(se_iN),
                 label = paste("ncol of assay", a))
    expect_equal(nrow(m), nrow(se_iN),
                 label = paste("nrow of assay", a))
  }
})

test_that("rowData contains V_Clamp for step-wise protocols", {
  data(se_iN)
  expect_true("V_Clamp" %in% colnames(rowData(se_iN)))
})

test_that("se_hAG has multiple assays", {
  data(se_hAG)
  expect_gte(length(assayNames(se_hAG)), 2L)
})

test_that("Well values follow plate-well format (letter + two digits)", {
  data(se_iN)
  wells <- colData(se_iN)$Well
  expect_true(all(grepl("^[A-P][0-9]{2}$", wells)))
})

test_that("no all-NA assay rows or columns", {
  data(se_iN)
  for (a in assayNames(se_iN)) {
    m <- assay(se_iN, a)
    expect_false(any(apply(m, 1, function(r) all(is.na(r)))),
                 label = paste("all-NA row in assay", a))
    expect_false(any(apply(m, 2, function(r) all(is.na(r)))),
                 label = paste("all-NA col in assay", a))
  }
})
