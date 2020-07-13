runSKATtraw <- function(phenodata, covariates, raw_file_path, window_size, window_shift, output_dir = "/scratch2/NSF_GWAS/Results/SKAT/", ncore="AllCores", SKAT_O = "OFF",
                        job_id, desired_sig_figs = 2, min_leading_0s = 2, upper_bound_leading_0s_null_model_par = 5, chunk_size,
                        switch_point = 4, plot = TRUE){ # get rid of switchpoint parameter once making function determine_switch_point?

  if(ncore=="AllCores") {
    ncore <- availableCores()
  }

  whole_genome_start <- proc.time()

  # Initial SKAT round w/o resampling, parallelized over SNP windows --------


  # null_model_noresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C")

  pre_allocated_SNP_windows <- pre_allocate(raw_file_path = raw_file_path,
                                            this_scaff_subset = this_scaff_subset,
                                            #window_list = window_list,
                                            window_size = window_size,
                                            window_shift = window_shift,
                                            pre_allocated_dir = "/scratch2/NSF_GWAS/SKAT_SLURMS/SKAT_MC_parallel/pre_allocated")

  #browser()

  master_output <- mtskat(this_phenotype = this_phenotype,
                          covariates = covariates,
                          raw_file_path = pre_allocated_SNP_windows[[1]][[3]],
                          this_scaff_subset = this_scaff_subset,
                          window_size = window_size,
                          window_shift = window_shift,
                          pre_allocated_SNP_windows = pre_allocated_SNP_windows,
                          chunk_size = chunk_size,
                          ncore = ncore)

  # Secondary SKAT round w/ resampling, parallelized over SNP windows -------


  max_leading_0s <- ceiling(-log(min(na.omit(master_output$`SKAT_p-val`)), base = 10)) - 1
  message(paste0("\n\nMax # leading 0s is: ", max_leading_0s, "\n\n"))

  for(leading_0s in min_leading_0s:min(switch_point, max_leading_0s)){


    add_to_master_output <- mtmcskat_SNPs(leading_0s = leading_0s,
                                          desired_sig_figs = desired_sig_figs,
                                          master_output = master_output,
                                          this_phenotype = this_phenotype,
                                          covariates = covariates,
                                          window_size = window_size,
                                          raw_file_path = pre_allocated_SNP_windows[[1]][[3]],
                                          ncore = ncore,
                                          pre_allocated_SNP_windows = pre_allocated_SNP_windows)

    if(add_to_master_output == "break") break

    # print(head(master_output))
    #
    # print(add_to_master_output)

    master_output <- rbind(master_output[-which(master_output$position %in% add_to_master_output$position), ],
                           add_to_master_output)
  }

  # Tertiary SKAT round w/ resampling, parallelized over null models --------

  lower_bound_leading_0s_null_model_par <- leading_0s

  for(leading_0s in lower_bound_leading_0s_null_model_par:upper_bound_leading_0s_null_model_par){

    # print(paste0("Running for", leading_0s, " leading zeroes"))

    if( leading_0s == upper_bound_leading_0s_null_model_par ){
      terminal_resampling <- TRUE
    }
    if( leading_0s < upper_bound_leading_0s_null_model_par ){
      terminal_resampling <- FALSE
    }
    if( leading_0s > upper_bound_leading_0s_null_model_par ){
      stop("Something is wrong. Attempting to run over range beyond terminal resampling.")
    }

    if(terminal_resampling == TRUE){
      lower_bound_p_val_for_MC_subset <- 0
    }
    if(terminal_resampling == FALSE){
      lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
    }
    upper_bound_p_val_for_MC_subset <- 0.1^leading_0s

    # How many permutations are needed for the desired # significant figures?
    n_permutations <- 10^(leading_0s+(desired_sig_figs))

    #browser()

    # How many jobs do we want to break this down into?
    # Makes sense to have one for each core to minimize communication
    n_jobs <- ncore

    #browser()

    #large_n_tests <- ceiling(n_permutations / n_jobs)

    max_permutations_per_job <- calculate_max_perm_per_core(
      nm_RAM_per_perm = benchmark_nm(phenotype = this_phenotype,
                                     covariates = covariates),
      RAM = size_RAM_partition(partition_factor = 2),
      ncore = ncore)

    #browser()

    # The number of permutations per null model should be such that
    # there is only one round of communication... unless we don't have
    # enough RAM for that, in which case we can divide jobs to give
    # each core twice as many (or more if needed)

    n_jobs <- ncore
    n_permutations_per_job <- ceiling(n_permutations / n_jobs)

    if(n_permutations_per_job >= max_permutations_per_job ){
      message("Memory constraint recognized")
      n_permutations_per_job <- max_permutations_per_job
      n_jobs <- ceiling(n_permutations / n_permutations_per_job)
    }
    #n_small_null_models <- n_permutations

    actual_n_permutations <- n_permutations_per_job * n_jobs

    #window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position

    window_list <- extract_windows_range_p(data = master_output,
                                           upper_bound = upper_bound_p_val_for_MC_subset,
                                           lower_bound = lower_bound_p_val_for_MC_subset)

    # Will be able to get rid of this line if we are quicker to clean out undesired rows withotu resampling
    window_list <- unique(window_list)

    message(paste0(Sys.time(), " - Running resampling for ", length(window_list),
                   " SNP windows with p-values between ",
                   upper_bound_p_val_for_MC_subset,
                   " and ", lower_bound_p_val_for_MC_subset,
                   " (", leading_0s, " leading zeroes) with at least ",
                   n_permutations, " permutations so empirical p-values",
                   " can be calculated to ", desired_sig_figs,
                   " significant figures..."))

    message(paste0("To run ", n_jobs, "jobs, each with ",
                   n_permutations_per_job, " permutations ",
                   "for a total of ", actual_n_permutations, " permutations\n"))

    if(length(window_list) == 0){
      message(paste0(Sys.time(),
                     " - No SNPs with p-vals within the range bounded from ",
                     lower_bound_p_val_for_MC_subset, " to ", upper_bound_p_val_for_MC_subset))
      #browser()
      next
    }

    new_pre_allocated_SNP_windows <- subset_list_of_lists(list_of_lists = pre_allocated_SNP_windows, desired_list = window_list, subindex = 1)
    #browser()
    timer <- proc.time()
    p_null_tallies <- future.apply::future_lapply(X = 1:n_jobs,
                                                  FUN = map_SKAT_nm,
                                                  pos_and_SNP_list = new_pre_allocated_SNP_windows,
                                                  window_size = window_size,
                                                  raw_file_path = new_pre_allocated_SNP_windows[[1]][[3]],
                                                  n_permutations = n_permutations_per_job,
                                                  resampling=TRUE)

    p_null_tallies <- dplyr::bind_rows(p_null_tallies)

    message(paste0(Sys.time(), " - Finished resampling up to ", actual_n_permutations, " permutations in",
                   (proc.time() - timer)[3], "s\n\n"))

    p_empirical_table <- p_empirical_from_tally(p_null_tallies = p_null_tallies,
                                                scaffold_ID = basename(raw_file_path))

    #add_to_master_output <- post_process_master_output(master_output = add_to_master_output[-1,])

    #browser()
    #add_to_master_output$`SKAT_p-val` <- as.numeric(as.character(output_resampled$`SKAT_p-val`))
    # print(dim(master_output))
    # print(master_output[1:6, 1:6])
    # print(dim(output_resampled))
    # print(output_resampled[1:6, 1:6])

    #browser()

    master_output <- rbind(master_output[-which(master_output$position %in% p_empirical_table$position), ],
                           p_empirical_table)
  }

  #browser()

  output_dir <- paste0(output_dir, "/", job_id)
  output_basename <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), ".csv")

  raw_out_name <- paste0(output_basename, ".csv")
  plot_out_name <- paste0(output_basename, ".png")

  message(paste0("Writing ", raw_out_name))

  fwrite(master_output, raw_out_name)

  if(plot==TRUE){

    compare_plots(master_output,
                  plot_out_name = plot_out_name,
                  scaffold_ID = pre_allocated_SNP_windows[[1]][[3]])

  }

  message(paste0("Done running in...",
                 print(proc.time() - whole_genome_start)[3],
                 "s"))

}
