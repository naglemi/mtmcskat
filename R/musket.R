mtskat <- function(this_phenotype,
                   covariates,
                   raw_file_path,
                   this_scaff_subset,
                   pre_allocated_SNP_windows,
                   pre_allocated_dir,
                   window_size,
                   window_shift,
                   chunk_size,
                   ncore){

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

  master_output <- loopSKAT(
    pre_allocated_SNP_windows = pre_allocated_SNP_windows,
    #grab_window = grab_window,
    window_size = window_size,
    window_shift = window_shift,
    #window_list = window_list,
    raw_file_path = raw_file_path,
    resampling = FALSE,
    null_model = null_model_noresample,
    n_permutations = NA,
    #this_scaff_subset = this_scaff_subset,
    backend = "furrr",
    ncore = ncore)

  master_output <- dplyr::bind_rows(master_output)
  master_output <- post_process_master_output(
    master_output = master_output[-1,])

  master_output
}
