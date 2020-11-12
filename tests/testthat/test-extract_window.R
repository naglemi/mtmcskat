context("SNP windows to be kernelized in SKAT")

data("sample_genodata")
data("sample_SNP_window")

test_that("SNP window is centered around desired position", {
  expect_equal(
    extract_window(this_position = 14460000,
                   window_size = 3000,
                   genodata = sample_genodata)$Position,
    14460000)
  }
)

test_that("SNP window is labeled with chromosome of origin", {
  expect_equal(
    extract_window(this_position = 14460000,
                   window_size = 3000,
                   genodata = sample_genodata)$Chr,
    10)
  }
)

test_that("Extracted SNP window is exactly the same as expected", {
  expect_equal(
    extract_window(this_position = 14460000,
                   window_size = 3000,
                   genodata = sample_genodata)$Z,
    sample_SNP_window$Z)
  }
)


test_that("Extracted SNP window and metadata is exactly the same as expected", {
  expect_equal(
    extract_window(this_position = 14460000,
                   window_size = 3000,
                   genodata = sample_genodata),
    sample_SNP_window)
  }
)

data("sample_pre_allocated_SNP_windows")

raw_file_path <- system.file("extdata",
                             "poplar_SNPs_Chr10_14460to14550kb.traw",
                             package = "SKATMCMT")

test_that("Length of pre-allocated window list is the same as expected", {
  expect_equal(
    length(pre_allocate(pre_allocated_dir = tempdir(),
                        raw_file_path = raw_file_path,
                        window_size = 3000,
                        window_shift = 1000)),
    90)
  }
)

test_that("Pre-allocated window list is exactly the same as expected", {
  expect_equal(
    pre_allocate(pre_allocated_dir = tempdir(),
                 raw_file_path = raw_file_path,
                 window_size = 3000,
                 window_shift = 1000),
    sample_pre_allocated_SNP_windows)
  }
)
