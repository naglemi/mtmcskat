mtmcskat_SNPs <- function(pre_allocated_SNP_windows,
                          this_phenotype,
                          covariates,
                          raw_file_path,
                          ncore){

  null_model <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates),
    n.Resampling = n_permutations,
    type.Resampling = "bootstrap")

  time_to_run_mapping <- proc.time()

  add_to_master_output <- future.apply::future_lapply(
    X = pre_allocated_SNP_windows,
    FUN = mappable_SKAT,
    #.options = future_options(packages = "SKAT"),
    #.progress = TRUE,
    #.cleanup = TRUE,
    #this_position = .x,
    #window_size = window_size,
    #this_scaff_subset = this_scaff_subset,
    raw_file_path = raw_file_path,
    null_model = null_model,
    resampling = TRUE,
    n_permutations = n_permutations,
    chunk = FALSE)

  message(paste0("Finished parallel run in ",
                 (proc.time() - time_to_run_mapping)[3],
                 "s"))

  add_to_master_output <- dplyr::bind_rows(add_to_master_output)
  add_to_master_output <- post_process_master_output(
    master_output = add_to_master_output)

  add_to_master_output
}
