context("Initial multi-threaded SKAT without resampling")

data("sample_mtskat_results")
data("sample_phenotype")
data("sample_covariates")
data("sample_pre_allocated_SNP_windows")

library(future)

old <- options(stringsAsFactors = FALSE)
options(mc.cores=2)

test_that("After pre-allocation, SNP windows are divided in desired n chunks", {
  expect_equal(length(chunk_windows(
    pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows,
    n_thread = 48)),
    48)
})

test_that("Multi-threaded SKAT produces output matching expected results", {
  expect_equal(mtskat(
    this_phenotype = sample_phenotype,
    covariates = sample_covariates,
    raw_file_path = sample_pre_allocated_SNP_windows[[1]][[3]],
    window_size = 3000,
    window_shift = 1000,
    pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows,
    n_thread = 2),
    sample_mtskat_results,
    tolerance = 1e-7) # Due to inequality that appears to be floating point error
})

on.exit(options(old), add = TRUE)
