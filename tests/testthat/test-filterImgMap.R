library(testthat)
library(ephacRTools)

# Helper: build a named map the same way localImgBrowserServer produces it
# Keys: "subdir|well.channel.class"  Values: URL strings
make_map <- function(...) {
  args <- list(...)
  setNames(unlist(args), names(args))
}

# ── filterImgMap — basic filtering ───────────────────────────────────────────

test_that("filterImgMap returns well->url for matching channel and class", {
  map <- make_map(
    "|A01.BF.class1"  = "blob:A01-BF-1",
    "|B02.BF.class1"  = "blob:B02-BF-1",
    "|A01.GFP.class1" = "blob:A01-GFP-1"
  )
  out <- filterImgMap(map, subdir = "", channel = "BF", img_class = "class1")
  expect_equal(sort(names(out)), c("A01", "B02"))
  expect_equal(out[["A01"]], "blob:A01-BF-1")
})

test_that("filterImgMap filters by img_class", {
  map <- make_map(
    "|A01.BF.class1" = "blob:A01-BF-1",
    "|A01.BF.class2" = "blob:A01-BF-2"
  )
  out1 <- filterImgMap(map, subdir = "", channel = "BF", img_class = "class1")
  out2 <- filterImgMap(map, subdir = "", channel = "BF", img_class = "class2")
  expect_equal(out1[["A01"]], "blob:A01-BF-1")
  expect_equal(out2[["A01"]], "blob:A01-BF-2")
})

test_that("filterImgMap filters by subdir (plate)", {
  map <- make_map(
    "plate1|A01.BF.class1" = "blob:p1-A01",
    "plate2|A01.BF.class1" = "blob:p2-A01"
  )
  out1 <- filterImgMap(map, subdir = "plate1", channel = "BF", img_class = "class1")
  out2 <- filterImgMap(map, subdir = "plate2", channel = "BF", img_class = "class1")
  expect_equal(out1[["A01"]], "blob:p1-A01")
  expect_equal(out2[["A01"]], "blob:p2-A01")
  expect_equal(length(out1), 1L)
})

test_that("filterImgMap returns empty vector when nothing matches", {
  map <- make_map("|A01.BF.class1" = "blob:url")
  out <- filterImgMap(map, subdir = "", channel = "GFP", img_class = "class1")
  expect_equal(length(out), 0L)
  expect_true(is.character(out))
})

test_that("filterImgMap keeps only first entry per well (dedup)", {
  # Two sites for A01, same key — first one wins
  map <- make_map(
    "|A01.BF.class1" = "blob:site1",
    "|A01.BF.class1" = "blob:site2"   # duplicate key: R keeps first
  )
  out <- filterImgMap(map, subdir = "", channel = "BF", img_class = "class1")
  expect_equal(length(out), 1L)
  expect_equal(names(out), "A01")
})

# ── filterImgMap — NULL / empty inputs ──────────────────────────────────────

test_that("filterImgMap returns empty on NULL map", {
  out <- filterImgMap(NULL, subdir = "")
  expect_equal(length(out), 0L)
})

test_that("filterImgMap returns empty on zero-length map", {
  out <- filterImgMap(character(0), subdir = "")
  expect_equal(length(out), 0L)
})

test_that("filterImgMap with NULL channel returns all channels for subdir", {
  map <- make_map(
    "|A01.BF.class1"  = "blob:BF",
    "|A01.GFP.class1" = "blob:GFP"
  )
  out <- filterImgMap(map, subdir = "", channel = NULL, img_class = "class1")
  # Both channels included, first per well wins
  expect_equal(length(out), 1L)   # deduped to one "A01"
})

test_that("filterImgMap with NULL img_class returns all classes", {
  map <- make_map(
    "|A01.BF.class1" = "blob:c1",
    "|B02.BF.class2" = "blob:c2"
  )
  out <- filterImgMap(map, subdir = "", channel = "BF", img_class = NULL)
  expect_equal(sort(names(out)), c("A01", "B02"))
})
