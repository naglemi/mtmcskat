# DEPRECATED

loopSKAT <- function(pre_allocated_SNP_windows, grab_window, window_size = NA, window_shift = NA, window_list = NA,
                     raw_file_path = NA, resampling = FALSE, null_model, n_permutations, parallelize_over = "windows", n_small_null_models = NA,
                     this_scaff_subset, ncore, backend="foreach", RAM_GB = 4, chunk = TRUE){

  if(parallelize_over=="windows" & backend=="furrr"){

    time_to_run_mapping <- proc.time()

    master_output <- future.apply::future_lapply(X = pre_allocated_SNP_windows,
                                                 FUN = mappable_SKAT,
                                                 #.options = future_options(packages = "SKAT"),
                                                 #.progress = TRUE,
                                                 #.cleanup = TRUE,
                                                 #this_position = .x,
                                                 #window_size = window_size,
                                                 #this_scaff_subset = this_scaff_subset,
                                                 raw_file_path = raw_file_path,
                                                 null_model = null_model,
                                                 resampling = resampling,
                                                 n_permutations = n_permutations,
                                                 chunk = chunk)
    message(paste0("Finished parallel run in ",
                 (proc.time() - time_to_run_mapping)[3],
                 "s"))

    master_output_df <- dplyr::bind_rows(master_output)

  }

  master_output
}
