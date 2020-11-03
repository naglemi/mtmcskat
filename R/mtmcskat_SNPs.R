mtmcskat_SNPs <- function(leading_0s,
                          master_output,
                          this_phenotype,
                          covariates,
                          window_size,
                          raw_file_path,
                          desired_sig_figs,
                          ncore = availableCores(),
                          pre_allocated_SNP_windows){
  upper_bound_p_val_for_MC_subset <- 0.1^leading_0s
  lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
  n_permutations <- 10^( leading_0s + desired_sig_figs )

  window_list <- extract_windows_range_p(
    data = master_output,
    upper_bound = upper_bound_p_val_for_MC_subset,
    lower_bound = lower_bound_p_val_for_MC_subset)

  #window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position

  message(paste0(Sys.time(), " - Running resampling for ", length(window_list),
                 " SNP windows with p-values from ",
                 upper_bound_p_val_for_MC_subset, " and ",
                 lower_bound_p_val_for_MC_subset,
                 " (", leading_0s, " leading zeroes) with ",
                 n_permutations, " permutations."))

  if(length(window_list) == 0){

    message(paste0(
      "No SNPs with p-vals within the range bounded from ",
      lower_bound_p_val_for_MC_subset, " to ",
      upper_bound_p_val_for_MC_subset))

    next
  }

  if(length(window_list) < ncore){
    # SWITCH MODE TO PARALLELIZE OVER NULL MODELS INTEAD OF WINDOWS
    # lower_bound_leading_0s_null_model_par <- leading_0s # Moved outside of this loop
    message(paste0(
      "Automatically switching modes for windows with p-vals <= 10^-",
      leading_0s,
      " because there are fewer windows (",
      length(window_list),
      ") than there are cores in use (",
      ncore,
      ").\n"))
    return("break")
  }

  null_model <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates),
    n.Resampling=n_permutations,
    type.Resampling="bootstrap")

  new_pre_allocated_SNP_windows <- subset_list_of_lists(
    list_of_lists = pre_allocated_SNP_windows,
    desired_list = window_list,
    subindex = 1)

  print(paste0("Numer of SNP windows outside of subset_list_of_lists is ",
               length(new_pre_allocated_SNP_windows)))

  if(length(new_pre_allocated_SNP_windows) == 0){
    stop("Error: There is no SNP data for windows in the list provided.")
  }

  add_to_master_output <- loopSKAT(
    pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
    window_size = window_size,
    #window_shift = window_shift,
    raw_file_path = raw_file_path,
    resampling = TRUE,
    null_model = null_model,
    n_permutations = n_permutations,
    chunk = FALSE,
    #this_scaff_subset = this_scaff_subset,
    backend = "furrr")

  add_to_master_output <- dplyr::bind_rows(add_to_master_output)
  add_to_master_output <- post_process_master_output(
    master_output = add_to_master_output)

  add_to_master_output
}
