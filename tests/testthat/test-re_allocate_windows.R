context("Selection of SNP windows, for resampling, based on initial p-values")

data("sample_mtskat_results")
data("sample_pre_allocated_SNP_windows")
data("sample_re_allocated_SNP_windows")

test_that("SNP windows selected are those with p-values in range of interest", {
  expect_equal(
    select_windows_range_p(
      data = sample_mtskat_results,
      upper_bound = 0.01,
      lower_bound = 0.001),
    c(14502000, 14503000, 14504000, 14505000))
  }
)

test_that("SNP windows extracted for p-value range are exactly as expected", {
  expect_equal(
    subset_list_of_lists(
      list_of_lists = sample_pre_allocated_SNP_windows,
      desired_list = c(14502000, 14503000, 14504000, 14505000),
      subindex = 1),
    sample_re_allocated_SNP_windows)
  }
)

test_that(paste("SNP windows extracted for p-value range are exactly as",
                "expected when passed one level up from `subset_list_of_lists`",
                "to `re_allocate_windows`"), {
  expect_equal(
    re_allocate_windows(
      data = sample_mtskat_results,
      upper_bound = 0.01,
      lower_bound = 0.001,
      pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows),
    sample_re_allocated_SNP_windows)
  }
)
