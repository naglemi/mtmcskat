context("Selection of SNP windows, for resampling, based on initial p-values")

data("small_mtskat_results")
data("small_pre_allocated_windows")
data("sample_re_allocated_SNP_windows")

test_that("SNP windows selected are those with p-values in range of interest", {
  expect_equal(
    select_windows_range_p(
      x = small_mtskat_results,
      upper_bound = 0.01,
      lower_bound = 0.001),
    14515000)
  }
)

modified_sample_mtskat_results <- small_mtskat_results
modified_sample_mtskat_results$`SKAT_p-val_resampled`[1] <- 0.009

test_that("Extra SNP windows for which initial guess was wrong were found", {
  expect_equal(
    select_windows_range_p(
      x = modified_sample_mtskat_results,
      upper_bound = 0.01,
      lower_bound = 0.001),
    c(14515000, 14491000))
  }
)

test_that("SNP windows extracted for p-value range are exactly as expected", {
  expect_equal(
    subset_list_of_lists(
      list_of_lists = small_pre_allocated_windows,
      desired_list = c(14515000, 14491000),
      subindex = 1),
    sample_re_allocated_SNP_windows)
  }
)

test_that(paste("SNP windows extracted for p-value range are exactly as",
                "expected when passed one level up from `subset_list_of_lists`",
                "to `re_allocate_windows`"), {
  expect_equal(
    re_allocate_windows(
      x = modified_sample_mtskat_results,
      upper_bound = 0.01,
      lower_bound = 0.001,
      pre_allocated_SNP_windows = small_pre_allocated_windows),
    sample_re_allocated_SNP_windows)
  }
)
