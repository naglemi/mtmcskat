context("SNP windows to be kernelized in SKAT")

small_genodata_path <- system.file("extdata",
                                   "poplar_200genotypes_14490to14520kb.traw",
                                   package = "mtmcskat")
small_genodata <- data.table::fread(small_genodata_path)

raw_file_path <- system.file("extdata",
                             "poplar_200genotypes_14490to14520kb.traw",
                             package = "mtmcskat")

temp <- paste0(tempdir(), "/mtmcskat_", stringi::stri_rand_strings(1, 5))

# Tests with unfiltered/unimputed SNP data --------------------------------
data("sample_SNP_window")

test_that("SNP window (unfiltered) is centered around desired position", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = FALSE,
                   remove_novar_SNPs = FALSE,
                   missing_cutoff = 0)$Position,
    145e5)
  }
)

test_that("SNP window (unfiltered) is labeled with chromosome of origin", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = FALSE,
                   remove_novar_SNPs = FALSE,
                   missing_cutoff = 0)$Chr,
    10)
  }
)

test_that("Extracted SNP window (unfiltered) is exactly the same as expected", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = FALSE,
                   remove_novar_SNPs = FALSE,
                   missing_cutoff = 0)$Z,
    sample_SNP_window$Z)
  }
)


test_that(paste("Extracted SNP window  (unfiltered) and metadata is",
                "exactly the same as expected"), {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = FALSE,
                   remove_novar_SNPs = FALSE,
                   missing_cutoff = 0),
    sample_SNP_window)
  }
)

test_that(paste("Length of pre-allocated window list (unfiltered)",
                "is the same as expected"), {
  expect_equal(
    length(pre_allocate(pre_allocated_dir = temp,
                        raw_file_path = raw_file_path,
                        window_size = 3000,
                        window_shift = 1000,
                        impute_to_mean = FALSE,
                        remove_novar_SNPs = FALSE,
                        missing_cutoff = 0)),
    30)
  }
)

data("small_pre_allocated_windows")

test_that(paste("Pre-allocated window list (unfiltered)",
                "is exactly the same as expected"), {
  expect_equal(
    pre_allocate(pre_allocated_dir = temp,
                 raw_file_path = raw_file_path,
                 window_size = 3000,
                 window_shift = 1000,
                 impute_to_mean = FALSE,
                 remove_novar_SNPs = FALSE,
                 missing_cutoff = 0),
    small_pre_allocated_windows)
  }
)




# Tests with filtered, imputed SNP data -----------------------------------
data("sample_SNP_window_processed")
temp <- paste0(tempdir(), "/mtmcskat_", stringi::stri_rand_strings(1, 5))

test_that("SNP window (filtered) is centered around desired position", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = TRUE,
                   remove_novar_SNPs = TRUE,
                   missing_cutoff = 0.15)$Position,
    145e5)
}
)

test_that("SNP window (filtered) is labeled with chromosome of origin", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = TRUE,
                   remove_novar_SNPs = TRUE,
                   missing_cutoff = 0.15)$Chr,
    10)
}
)

test_that("Extracted SNP window (filtered) is exactly the same as expected", {
  expect_equal(
    extract_window(this_position = 145e5,
                   window_size = 3000,
                   genodata = small_genodata,
                   impute_to_mean = TRUE,
                   remove_novar_SNPs = TRUE,
                   missing_cutoff = 0.15)$Z,
    sample_SNP_window_processed$Z)
}
)


test_that(paste("Extracted SNP window  (filtered) and metadata is",
                "exactly the same as expected"), {
                  expect_equal(
                    extract_window(this_position = 145e5,
                                   window_size = 3000,
                                   genodata = small_genodata,
                                   impute_to_mean = TRUE,
                                   remove_novar_SNPs = TRUE,
                                   missing_cutoff = 0.15),
                    sample_SNP_window_processed)
                }
)

test_that(paste("Length of pre-allocated window list (unfiltered)",
                "is the same as expected"), {
                  expect_equal(
                    length(pre_allocate(pre_allocated_dir = temp,
                                        raw_file_path = raw_file_path,
                                        window_size = 3000,
                                        window_shift = 1000,
                                        impute_to_mean = TRUE,
                                        remove_novar_SNPs = TRUE,
                                        missing_cutoff = 0.15)),
                    30)
                }
)

data("small_pre_allocated_windows_processed")

test_that(paste("Pre-allocated window list (unfiltered)",
                "is exactly the same as expected"), {
                  expect_equal(
                    pre_allocate(pre_allocated_dir = temp,
                                 raw_file_path = raw_file_path,
                                 window_size = 3000,
                                 window_shift = 1000,
                                 impute_to_mean = TRUE,
                                 remove_novar_SNPs = TRUE,
                                 missing_cutoff = 0.15),
                    small_pre_allocated_windows_processed)
                }
)


