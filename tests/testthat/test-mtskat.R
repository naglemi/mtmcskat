context("Initial multi-threaded SKAT without resampling")

data("small_mtskat_results")
data("small_phenotype")
data("small_covariates")
data("small_pre_allocated_windows")

library(future)

old <- options(stringsAsFactors = FALSE)
options(mc.cores=2)

set.seed(5)

test_that("After pre-allocation, SNP windows are divided in desired n chunks", {
  expect_equal(length(chunk_windows(
    pre_allocated_SNP_windows = small_pre_allocated_windows,
    n_thread = 3)),
    3)
})

test_that("Multi-threaded SKAT produces output matching expected results", {
  expect_equal(mtskat(
    this_phenotype = small_phenodata,
    covariates = small_covariates,
    raw_file_path = small_pre_allocated_windows[[1]][[3]],
    pre_allocated_SNP_windows = small_pre_allocated_windows,
    n_thread = 2),
    small_mtskat_results,
    tolerance = 1e-7) # Due to inequality that appears to be floating point error
})

on.exit(options(old), add = TRUE)
