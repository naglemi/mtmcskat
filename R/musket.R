chunk_windows <- function(pre_allocated_SNP_windows,
                          n_core){

  chunk_size <- length(pre_allocated_SNP_windows) / n_core

  message(paste0(Sys.time(), " - Chunking list of length ",
                 length(pre_allocated_SNP_windows),
                 " into blocks each with no more than",
                 chunk_size,
                 " SNP windows each"))

  pre_allocated_SNP_windows <- split(
    pre_allocated_SNP_windows,
    (seq_along(pre_allocated_SNP_windows) - 1) %/% chunk_size)

  message(paste0("...resulting in ",
                 length(pre_allocated_SNP_windows),
                 " blocks"))

  pre_allocated_SNP_windows
}

mtskat <- function(this_phenotype,
                   covariates,
                   raw_file_path,
                   this_scaff_subset,
                   pre_allocated_SNP_windows,
                   pre_allocated_dir,
                   window_size,
                   window_shift,
                   chunk_size,
                   n_core){

  null_model_noresample <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates), out_type="C")

  if(!exists("pre_allocated_SNP_windows")){

      pre_allocated_SNP_windows <- pre_allocate(
        raw_file_path = raw_file_path,
        this_scaff_subset = this_scaff_subset,
        #window_list = window_list,
        window_size = window_size,
        window_shift = window_shift,
        pre_allocated_dir = pre_allocated_dir)
  }

  pre_allocated_SNP_windows <- chunk_windows(
    pre_allocated_SNP_windows = pre_allocated_SNP_windows,
    n_core = n_core)

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
