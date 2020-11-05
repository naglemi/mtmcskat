context("Empirical p-values")

data("sample_covariates")
data("sample_null_model")
data("sample_pre_allocated_SNP_windows")
data("sample_phenotype")

test_that("Empirical p-val is equal to earlier result with same null model", {
  expect_equal(
    calculate_SKAT_empirical_p(
      Z = sample_pre_allocated_SNP_windows[[1]][[2]],
      n_permutations = 1000,
      null_model = sample_null_model,
      return_all_p = FALSE),
    0.01298701)
  }
)

test_that("Vector of empirical p-values is equal to # permutations", {
  expect_equal(
    length(
      calculate_SKAT_empirical_p(
        Z = sample_pre_allocated_SNP_windows[[1]][[2]],
        n_permutations = 1000,
        null_model = sample_null_model,
        return_all_p = TRUE)),
    1000)
  }
)

