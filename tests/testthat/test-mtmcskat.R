context("MTMCSKAT with multithreading over SNP windows and over null models")

data("small_mtskat_results")
data("small_phenodata")
data("small_covariates")
data("small_pre_allocated_windows")
data("sample_mtmcskat_results")
data("sample_mtmcskatnm1_results")
data("sample_mtmcskatnm2_results")

old <- options(stringsAsFactors = FALSE)
options(mc.cores=2)

set.seed(5)

test_that("MTMCSKAT threading on SNP windows gives expected results", {
  expect_equal(mtmcskat_SNPs(
    this_phenotype = small_phenodata,
    covariates = small_covariates,
    n_permutations = 500,
    pre_allocated_SNP_windows = small_pre_allocated_windows[2:4],
    scaffold_ID = small_pre_allocated_windows[[1]][[3]],
    n_thread = 2),
    sample_mtmcskat_results)
})

test_that("MTMCSKAT threading on models in 1 batch gives expected results", {
  expect_equal(mtmcskat_NullModels(
    this_phenotype = small_phenodata,
    covariates = small_covariates,
    n_permutations = 500,
    n_thread = 2,
    max_permutations_per_job = 251,
    pre_allocated_SNP_windows = small_pre_allocated_windows[2:4],
    scaffold_ID = small_pre_allocated_windows[[1]][[3]]),
    sample_mtmcskatnm1_results)
})

test_that("MTMCSKAT threading on models in 2 batches gives expected results", {
  expect_equal(mtmcskat_NullModels(
    this_phenotype = small_phenodata,
    covariates = small_covariates,
    n_permutations = 500,
    n_thread = 2,
    max_permutations_per_job = 249,
    pre_allocated_SNP_windows = small_pre_allocated_windows[2:4],
    scaffold_ID = small_pre_allocated_windows[[1]][[3]]),
    sample_mtmcskatnm2_results)
})

on.exit(options(old), add = TRUE)
