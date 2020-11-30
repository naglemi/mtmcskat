context("Empirical p-values")

data("small_covariates")
data("small_pre_allocated_windows")
data("small_phenodata")

set.seed(5)

sample_null_model <- SKAT::SKAT_Null_Model(
  small_phenodata ~ 1 + as.matrix(small_covariates), out_type="C",
  n.Resampling = 1000)

test_that("Empirical p-val is equal to earlier result with same null model", {
  expect_equal(
    calculate_SKAT_empirical_p(
      Z = small_pre_allocated_windows[[18]][[2]],
      n_permutations = 1000,
      null_model = sample_null_model,
      return_all_p = FALSE),
    0.019)
  }
)

test_that("Vector of empirical p-values is equal to # permutations", {
  expect_equal(
    length(
      calculate_SKAT_empirical_p(
        Z = small_pre_allocated_windows[[18]][[2]],
        n_permutations = 1000,
        null_model = sample_null_model,
        return_all_p = TRUE)),
    1000)
  }
)

test_that("n_permutation must match # permutations in null model submitted", {
  expect_error(
    calculate_SKAT_empirical_p(
      Z = sample_pre_allocated_windows[[31]][[2]],
      n_permutations = 500,
      null_model = sample_null_model,
      return_all_p = TRUE))
  }
)

test_that(paste("Empirical p-val passed up one level to SKAT_one_window",
                "is as expected"), {
                  expect_equal(
                    SKAT_one_window(
                      this_position = small_pre_allocated_windows[[18]][[1]],
                      Z = small_pre_allocated_windows[[18]][[2]],
                      scaffold_ID = small_pre_allocated_windows[[18]][[3]],
                      n_permutations = 1000,
                      null_model = sample_null_model,
                      resampling = TRUE,
                      return_all_p = FALSE)[4],
                    0.019)
                }
)

test_that(paste("Empirical p-val passed up two levels to mappable_SKAT",
                "is as expected"), {
                  expect_equal(
                    mappable_SKAT(
                      pos_and_SNPs = small_pre_allocated_windows[[18]],
                      scaffold_ID = small_pre_allocated_windows[[18]][[3]],
                      null_model = sample_null_model,
                      resampling = TRUE,
                      n_permutations = 1000,
                      chunk = FALSE)$V4,
                    0.019)
                }
)

