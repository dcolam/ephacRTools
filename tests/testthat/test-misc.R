library(testthat)
library(ephacRTools)

# ── dround ──────────────────────────────────────────────────────────────────

test_that("dround returns numeric output", {
  expect_true(is.numeric(dround(0.00002345)))
  expect_true(is.numeric(dround(c(0.1, 100, 0.0001))))
})

test_that("dround trims small numbers to correct significant digits", {
  expect_equal(dround(0.00002345, digits = 3), 0.0000234)
  expect_equal(dround(0.001234,   digits = 2), 0.0012)
})

test_that("dround handles numbers greater than 1", {
  expect_equal(dround(554356, digits = 3), 554000)
  expect_equal(dround(12.56,  digits = 3), 12.6)
})

test_that("dround preserves vector length", {
  x <- c(0.00002345, 554356, 12.56)
  expect_equal(length(dround(x)), length(x))
})

test_that("dround works on data.frames (numeric columns only)", {
  df <- data.frame(a = c(0.001234, 0.005678), b = c("x", "y"),
                   stringsAsFactors = FALSE)
  out <- dround(df, digits = 2)
  expect_equal(out$a, c(0.0012, 0.0057))
  expect_equal(out$b, c("x", "y"))   # character column unchanged
})

test_that("dround with roundGreaterThan1 = FALSE leaves large numbers alone", {
  expect_equal(dround(554356, digits = 3, roundGreaterThan1 = FALSE), 554356)
  expect_equal(dround(0.001234, digits = 2, roundGreaterThan1 = FALSE), 0.0012)
})


# ── grepGene ─────────────────────────────────────────────────────────────────

test_that("grepGene finds genes by symbol in composite names", {
  x <- c("ENSG00000000003.TSPAN6", "ENSG00000000005.TNMD",
         "ENSG00000000419.DPM1")
  expect_equal(grepGene(x, "TSPAN6"), "ENSG00000000003.TSPAN6")
  expect_equal(sort(grepGene(x, c("TSPAN6", "TNMD"))),
               sort(c("ENSG00000000003.TSPAN6", "ENSG00000000005.TNMD")))
})

test_that("grepGene returns character(0) for no match", {
  x <- c("ENSG00000000003.TSPAN6", "ENSG00000000005.TNMD")
  expect_equal(grepGene(x, "NOTEXIST"), character(0))
})

test_that("grepGene returns exact match when gene is already in vector", {
  x <- c("TSPAN6", "TNMD", "DPM1")
  expect_equal(grepGene(x, "TSPAN6"), "TSPAN6")
})

test_that("grepGene is case-insensitive by default", {
  x <- c("ENSG00000000003.TSPAN6", "ENSG00000000005.TNMD")
  expect_equal(grepGene(x, "tspan6"), "ENSG00000000003.TSPAN6")
})

test_that("grepGene subsets SummarizedExperiment by row", {
  data(se_iN)
  result <- grepGene(se_iN, rownames(se_iN)[1:2])
  expect_s4_class(result, "SingleCellExperiment")
  expect_equal(nrow(result), 2L)
})


# ── %||% ─────────────────────────────────────────────────────────────────────

test_that("%||% returns right-hand side when left is NULL", {
  expect_equal(NULL %||% "fallback", "fallback")
  expect_equal(NULL %||% 42L,        42L)
  expect_equal(NULL %||% NULL,       NULL)
})

test_that("%||% returns left-hand side when it is not NULL", {
  expect_equal("value" %||% "fallback", "value")
  expect_equal(0       %||% 99,          0)
  expect_equal(FALSE   %||% TRUE,        FALSE)
})
