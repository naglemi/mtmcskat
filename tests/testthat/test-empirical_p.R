context("Empirical p-values")

data("sample_covariates")
data("sample_null_model")
data("sample_pre_allocated_SNP_windows")
data("sample_phenotype")

test_that("Empirical p-val is equal to earlier result with same null model", {
  expect_equal(
    calculate_SKAT_empirical_p(
      Z = sample_pre_allocated_SNP_windows[[31]][[2]],
      n_permutations = 1000,
      null_model = sample_null_model,
      return_all_p = FALSE),
    0.1008991)
  }
)

test_that("Vector of empirical p-values is equal to # permutations", {
  expect_equal(
    length(
      calculate_SKAT_empirical_p(
        Z = sample_pre_allocated_SNP_windows[[31]][[2]],
        n_permutations = 1000,
        null_model = sample_null_model,
        return_all_p = TRUE)),
    1000)
  }
)

test_that("n_permutation must match # permutations in null model submitted", {
  expect_error(
    calculate_SKAT_empirical_p(
      Z = sample_pre_allocated_SNP_windows[[31]][[2]],
      n_permutations = 500,
      null_model = sample_null_model,
      return_all_p = TRUE))
  }
)

test_that(paste("Empirical p-val passed up one level to SKAT_one_window",
                "is as expected"), {
                  expect_equal(
                    SKAT_one_window(
                      this_position = sample_pre_allocated_SNP_windows[[31]][[1]],
                      Z = sample_pre_allocated_SNP_windows[[31]][[2]],
                      scaffold_ID = sample_pre_allocated_SNP_windows[[31]][[3]],
                      n_permutations = 1000,
                      null_model = sample_null_model,
                      resampling = TRUE,
                      return_all_p = FALSE)[4],
                    0.1008991)
                }
)

test_that(paste("Empirical p-val passed up two levels to mappable_SKAT",
                "is as expected"), {
                  expect_equal(
                    mappable_SKAT(
                      pos_and_SNPs = sample_pre_allocated_SNP_windows[[31]],
                      scaffold_ID = sample_pre_allocated_SNP_windows[[31]][[3]],
                      null_model = sample_null_model,
                      resampling = TRUE,
                      n_permutations = 1000,
                      chunk = FALSE)$V4,
                    0.1008991)
                }
)

