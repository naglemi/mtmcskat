#' Wrapper for SKAT sans resampling, with multi-threading over SNP windows
#'
#' @param this_phenotype Phenotype data that has been read into R from the
#'   standard PLINK phenotype file format
#'   (\href{https://www.cog-genomics.org/plink2/formats}{PLINK documentation})
#' @param covariates Covariates data that has been read into R from a
#'   comma-delimited text file with a header and one row for each genotype, in
#'   the same order as the phenotype file.
#' @param pre_allocated_SNP_windows A list of SNP windows along with their
#'   positions and scaffolds, prepared using \code{\link{pre_allocate}}
#' @inheritParams pre_allocate
#' @inheritParams calculate_SKAT_empirical_p
#' @inheritParams SKAT_one_window
#' @inheritParams chunk_windows
#'
#' @return A dataframe with columns for scaffold ID, SNP window position,
#'   p-values from the model used in SKAT without resampling, and an additional
#'   column for empirical p-values (all NA at this stage since resampling is not
#'   run in `mtskat`)
#' @export
#'
#' @examples
#' \dontrun{
#' data("sample_phenotype")
#' data("sample_covariates")
#' data("sample_pre_allocated_SNP_windows")
#'
#' mtskat(
#'   this_phenotype = sample_phenotype,
#'   covariates = sample_covariates,
#'   raw_file_path = sample_pre_allocated_SNP_windows[[1]][[3]],
#'   window_size = 3000,
#'   window_shift = 1000,
#'   pre_allocated_SNP_windows = sampel_pre_allocated_SNP_windows,
#'   n_core = 2)
#' }
mtskat <- function(this_phenotype,
                   covariates,
                   raw_file_path,
                   pre_allocated_SNP_windows,
                   pre_allocated_dir,
                   window_size,
                   window_shift,
                   n_thread){

  null_model_noresample <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates), out_type="C")

  if(!exists("pre_allocated_SNP_windows")){

      pre_allocated_SNP_windows <- pre_allocate(
        raw_file_path = raw_file_path,
        window_size = window_size,
        window_shift = window_shift,
        pre_allocated_dir = pre_allocated_dir)
  }

  message(paste("Memory taken up by pre-allocated SNP windows before chunking is",
                object.size(pre_allocated_SNP_windows)/1e6,
                "MB"))

  pre_allocated_SNP_windows <- chunk_windows(
    pre_allocated_SNP_windows = pre_allocated_SNP_windows,
    n_thread = n_thread)

  message(paste("Memory taken up by pre-allocated SNP windows AFTER chunking is",
                object.size(pre_allocated_SNP_windows)/1e6,
                "MB"))

  time_to_run_mapping <- proc.time()

  print("Size of all objects in environment, in MB is: ")
  print(sort( sapply(ls(),function(x){object.size(get(x))/1e6})))

  master_output <- future.apply::future_lapply(X = pre_allocated_SNP_windows,
                                               FUN = mappable_SKAT,
                                               scaffold_ID = raw_file_path,
                                               null_model = null_model_noresample,
                                               resampling = FALSE,
                                               n_permutations = NA,
                                               chunk = TRUE)
  message(paste0("Finished parallel run in ",
                 (proc.time() - time_to_run_mapping)[3],
                 "s"))

  master_output <- dplyr::bind_rows(master_output)
  master_output <- post_process_master_output(
    master_output = master_output[-1,])

  master_output
}
