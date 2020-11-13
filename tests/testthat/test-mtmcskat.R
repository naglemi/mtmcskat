
data("sample_mtskat_results")
data("sample_phenotype")
data("sample_covariates")
data("sample_pre_allocated_SNP_windows")
data("sample_mtmcskat_results")
data("sample_mtmcskatnm1_results")
data("sample_mtmcskatnm2_results")

old <- options(stringsAsFactors = FALSE)
options(mc.cores=2)

set.seed(5)

test_that("MTMCSKAT threading on SNP windows gives expected results", {
  expect_equal(mtmcskat_SNPs(
    this_phenotype = sample_phenotype,
    covariates = sample_covariates,
    n_permutations = 500,
    pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
    scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]]),
    sample_mtmcskat_results)
})

test_that("MTMCSKAT threading on models in 1 batch gives expected results", {
  expect_equal(mtmcskat_NullModels(
    this_phenotype = sample_phenotype,
    covariates = sample_covariates,
    n_permutations = 500,
    n_thread = 2,
    max_permutations_per_job = 251,
    pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
    scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]]),
    sample_mtmcskatnm1_results)
})

test_that("MTMCSKAT threading on models in 2 batches gives expected results", {
  expect_equal(mtmcskat_NullModels(
    this_phenotype = sample_phenotype,
    covariates = sample_covariates,
    n_permutations = 500,
    n_thread = 2,
    max_permutations_per_job = 249,
    pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
    scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]]),
    sample_mtmcskatnm2_results)
})

on.exit(options(old), add = TRUE)
