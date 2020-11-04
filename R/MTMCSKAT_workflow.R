#' Run an automated MTMC-SKAT job
#'
#' This function provides a one-liner to run MTMC-SKAT over a scaffold.
#'
#' For whole-genome analysis, this function (accesible by command line) is to be parallelized over jobs submitted in a batch query system on a high-performance cluster, or looped over scaffolds if running on a single machine.
#'
#' @param phenodata string for filepath for phenotype file, with labeled columns for FID, IID, and trait.
#' @param covariates string for filepath for covariate file, with ordered covariate values in each column and no header
#' @param raw_file_path string for filepath to .traw file (see [PLINK file format reference](https://www.cog-genomics.org/plink/2.0/formats))
#' @param window_size numeric, the size of each SNP window for which the SKAT kernalization is to be performed, in base pairs
#' @param window_shift numeric, the number of base pairs over which each rolling window will slide; in other terms, the distance between the start (or end) positions of adjacent overlapping windows
#' @param output_dir string for directory where results will be saved
#' @param n_core numeric, maximum number of cores over which parallelization is allowed
#' @param job_id string, identifier to label the job; this identifier will go into output filename
#' @param desired_sig_figs number of significant figures desired for
#' @param min_leading_0s numeric, threshold x to begin resampling for SNP windows with p-values below 10^-x; default value is 2
#' @param max_accuracy numeric, limit x to end resampling for SNP windows with p-values below 10^-x, usually due to computational cost
#' @param chunk_size DEPRECATED?
#' @param switch_point numeric, limit x at which SNP windows with p-values below 10^-x must be testing by parallelizing over null models rather than SNP windows (may be set by user due to limitations on RAM preventing production of large null models)
#' @param plot if TRUE, produce Manhattan plot of results
#'
#' @return None; outputs are saved to the user-specified output directory, with the user-specified ```job_id```
#' @export
#'
#' @examples
MTMCSKAT_workflow <- function(phenodata, covariates, raw_file_path, window_size,
                              window_shift, output_dir, pre_allocated_dir,
                              n_core="AllCores", job_id, desired_sig_figs = 2,
                              min_leading_0s = 2,
                              max_accuracy = 5,
                              chunk_size, switch_point = 4, plot = TRUE){ # get rid of switchpoint parameter once making function determine_switch_point?

  set_accuracy <- function(master_output,
                           max_accuracy){

    max_leading_0s <- ceiling(-log(
      min(na.omit(master_output$`SKAT_p-val`)),
      base = 10)) - 1

    message(paste0("\n\nMax # leading 0s is: ", max_leading_0s, "\n\n"))

    if (max_accuracy != "Auto"){
      message(paste0("\n\nUser-defined accuracy is 10^-", max_accuracy, "\n\n"))
      if (max_leading_0s < max_accuracy){
        message(paste0("This level of accuracy may need be needed since p-values",
                       "from mtskat have no more than ", max_leading_0s,
                       "leading zeros"))
      }
      if (max_leading_0s > max_accuracy){
        message(paste0("This accuracy may not be enough to accurately calculate ",
                       "all empirical p-values since p-values from mtskat have ",
                       "as many as ", max_leading_0s, "leading zeros"))
      }
    }
    if (max_accuracy == "Auto"){
      max_accuracy <- max_leading_0s
    }

    max_accuracy

  }

  determine_whether_final_resampling <- function(leading_0s,
                                                 max_accuracy){
    if( leading_0s == max_accuracy ){
      terminal_resampling <- TRUE
    }
    if( leading_0s < max_accuracy ){
      terminal_resampling <- FALSE
    }
    if( leading_0s > max_accuracy ){
      stop(paste0("Something is wrong. ",
                  "Attempting to run over range beyond terminal resampling."))
    }

    terminal_resampling

  }

  determine_subset_pval_boundaries <- function(leading_0s,
                                               desired_sig_figs,
                                               terminal_resampling){

    if(terminal_resampling == TRUE){
      lower_bound_p_val_for_MC_subset <- 0
    }
    if(terminal_resampling == FALSE){
      lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
    }
    upper_bound_p_val_for_MC_subset <- 0.1^leading_0s

    boundaries <- list("upper" = upper_bound_p_val_for_MC_subset,
                       "lower" = lower_bound_p_val_for_MC_subset)

    boundaries
  }

  determine_n_permutations <- function(leading_0s,
                                       desired_sig_figs){
    n_permutations <- 10^( leading_0s + desired_sig_figs )

    n_permutations
  }

  re_allocate_windows <- function(master_output,
                                  upper_bound_p_val_for_MC_subset,
                                  lower_bound_p_val_for_MC_subset,
                                  pre_allocated_SNP_windows,
                                  window_list){

    window_list <- extract_windows_range_p(
      data = master_output,
      upper_bound = upper_bound_p_val_for_MC_subset,
      lower_bound = lower_bound_p_val_for_MC_subset)

    if(length(window_list) == 0){

      message(paste0(Sys.time(),
                     " - No SNPs with p-vals within the range bounded from ",
                     lower_bound_p_val_for_MC_subset, " to ",
                     upper_bound_p_val_for_MC_subset))

      return("None")
    }

    new_pre_allocated_SNP_windows <- subset_list_of_lists(
      list_of_lists = pre_allocated_SNP_windows,
      desired_list = window_list,
      subindex = 1)

    print(paste0("Number of SNP windows outside of subset_list_of_lists is ",
                 length(new_pre_allocated_SNP_windows)))

    if(length(new_pre_allocated_SNP_windows) == 0){
      stop("Error: There is no SNP data for windows in the list provided.")
    }

    new_pre_allocated_SNP_windows

  }

  check_if_more_threads_than_windows <- function(
    pre_allocated_SNP_windows,
    n_core){
    # If there are fewer SNP windows than cores, we should parallelize over nms
    if(length(new_pre_allocated_SNP_windows) < n_core){
      # SWITCH MODE TO PARALLELIZE OVER NULL MODELS INTEAD OF WINDOWS

      message(paste0(
        "Automatically switching modes for windows with p-vals <= 10^-",
        leading_0s,
        " because there are fewer windows (",
        length(pre_allocated_SNP_windows),
        ") than there are cores in use (",
        n_core,
        ").\n"))
      return(TRUE)
    }

    if(length(new_pre_allocated_SNP_windows) >= n_core){
      return(FALSE)
    }

  }

  check_if_null_models_fit_in_RAM_per_thread <- function(n_permutations,
                                                       RAM_per_permutation,
                                                       RAM_per_thread){
    # to parallelize over SNP windows, we need to be able to fit a copy
    # of the (same) null model for each thread
    RAM_per_thread_needed_for_mtmcskat_SNPs <-
      RAM_per_permutation * RAM_per_thread

    # If the amount of RAM per core is less than this,
    # parallelize over null models instead of SNP windows.
    # Note in this latter case, null models for each core will be different!

    if(RAM_per_thread_needed_for_mtmcskat_SNPs < RAM_per_thread){
      null_models_fit_in_RAM_per_thread <- TRUE
    }

    if(RAM_per_thread_needed_for_mtmcskat_SNPs > RAM_per_thread){
      null_models_fit_in_RAM_per_thread <- FALSE
    }

    null_models_fit_in_RAM_per_thread

  }

  arrange_jobs_NullModel_multithreading <-
    function(n_core,
             n_permutations,
             max_permutations_per_job){
      # How many jobs do we want to break this down into?
      # Makes sense to have one for each core to minimize communication
      n_jobs <- n_core

      # The number of permutations per null model should be such that
      # there is only one round of communication... unless we don't have
      # enough RAM for that, in which case we can divide jobs to give
      # each core twice as many (or more if needed)

      n_total_permutations_per_thread <- ceiling(n_permutations / n_core)

      if(n_total_permutations_per_thread >= max_permutations_per_job ){
        message(paste0("Memory constraint recognized; ",
                       "Unable to multithread over all null models ",
                       "simultaneously."))

        n_permutations_per_job <- max_permutations_per_job
        n_jobs <- ceiling(n_permutations / n_permutations_per_job)
      }

      if(n_total_permutations_per_thread < max_permutations_per_job ){
        message(paste0("Memory is enough to multithread over all null ",
                       "models simultaneously."))

        n_permutations_per_job <- n_total_permutations_per_thread
        n_jobs <- n_core
      }

      # Due to rounding up and using same n permutations in null model for
      # each core, we may have more permutations than originally planned
      actual_n_permutations <- n_permutations_per_job * n_jobs

      job_details <- list(
        "n_permutations_per_job" = n_permutations_per_job,
        "n_jobs" = n_jobs,
        "actual_n_permutations" = actual_n_permutations)

      job_details
    }

  this_phenotype <- unlist(read.csv(phenodata,
                                    sep = "\t",
                                    colClasses = c("character",
                                                   "character",
                                                   "numeric"))[, 3])

  covariates <- read.csv(covariates,
                         sep = "\t",
                         header = FALSE)

  if(n_core=="AllCores") {
    n_core <- future::availableCores()
  }

  whole_genome_start <- proc.time()

  # Initial SKAT round w/o resampling, parallelized over SNP windows --------

  pre_allocated_SNP_windows <- pre_allocate(
    raw_file_path = raw_file_path,
    this_scaff_subset = this_scaff_subset,
    #window_list = window_list,
    window_size = window_size,
    window_shift = window_shift,
    pre_allocated_dir = pre_allocated_dir)

  master_output <- mtskat(
    this_phenotype = this_phenotype,
    covariates = covariates,
    raw_file_path = pre_allocated_SNP_windows[[1]][[3]],
    this_scaff_subset = this_scaff_subset,
    window_size = window_size,
    window_shift = window_shift,
    pre_allocated_SNP_windows = pre_allocated_SNP_windows,
    chunk_size = chunk_size,
    n_core = n_core)

  # Secondary SKAT round w/ resampling, multithread over SNP windows or NMs-----

  max_accuracy <- set_accuracy(master_output = master_output,
                               max_accuracy = max_accuracy)

  ## A null model with HOW MUCH resampling can fit into the RAM allotted
  ## for a single thread?

  RAM_per_permutation <- benchmark_nm(phenotype = this_phenotype,
                                      covariates = covariates,
                                      benchmark_size = 1000)

  # This parameter only used if multithreading over null models
  max_permutations_per_job <- calculate_max_perm_per_core(
    nm_RAM_per_perm = RAM_per_permutation,
    RAM = size_RAM_wiggle(wiggle_factor = 2),
    n_core = n_core)

  RAM_per_thread <- benchmarkme::get_ram() / n_core

  #for(leading_0s in min_leading_0s:min(switch_point, max_accuracy)){
  for(leading_0s in min_leading_0s:max_accuracy){

    n_permutations <- determine_n_permutations(
      leading_0s = leading_0s,
      desired_sig_figs = desired_sig_figs)

    terminal_resampling <- determine_whether_final_resampling(
      leading_0s = leading_0s,
      max_accuracy = max_accuracy)

    boundaries <- determine_subset_pval_boundaries(
      leading_0s = leading_0s,
      desired_sig_figs = desired_sig_figs,
      terminal_resampling = terminal_resampling)

    new_pre_allocated_SNP_windows <- re_allocate_windows(
      master_output = master_output,
      upper_bound_p_val_for_MC_subset = boundaries$upper,
      lower_bound_p_val_for_MC_subset = boundaries$lower,
      pre_allocated_SNP_windows = pre_allocated_SNP_windows,
      window_list = window_list)

    if(new_pre_allocated_SNP_windows == "None") next

    more_threads_than_windows <- check_if_more_threads_than_windows(
      pre_allocated_SNP_windows = pre_allocated_SNP_windows,
      n_core = n_core)

    null_models_fit_in_RAM_per_thread <-
      check_if_null_models_fit_in_RAM_per_thread(
        n_permutations = n_permutations,
        RAM_per_permutation = RAM_per_permutation,
        RAM_per_thread = RAM_per_thread)

    message(paste0(Sys.time(), " - Running resampling for ",
                   length(new_pre_allocated_SNP_windows),
                   " SNP windows with p-values between ",
                   boundaries$upper, " and ",
                   boundaries$lower,
                   " (", leading_0s, " leading zeroes) with at least ",
                   n_permutations, " permutations so empirical p-values",
                   " can be calculated to ", desired_sig_figs,
                   " significant figures..."))

    if(null_models_fit_in_RAM_per_thread & !more_threads_than_windows){
      message("Multithreading over SNP windows")

      add_to_master_output <- mtmcskat_SNPs(
        this_phenotype = this_phenotype,
        covariates = covariates,
        raw_file_path = pre_allocated_SNP_windows[[1]][[3]],
        n_core = n_core,
        pre_allocated_SNP_windows = pre_allocated_SNP_windows)
    }

    if(!null_models_fit_in_RAM_per_thread | more_threads_than_windows){
      message("Multithreading over null models")

      job_details <- arrange_jobs_NullModel_multithreading(
        n_core = n_core,
        n_permutations = n_permutations,
        max_permutations_per_job = max_permutations_per_job
      )

      message(paste0("To run ", job_details$n_jobs, "jobs, each with ",
                     job_details$n_permutations_per_job, " permutations ",
                     "for a total of ", job_details$actual_n_permutations,
                     " permutations\n"))

      timer <- proc.time()

      p_null_tallies <- future.apply::future_lapply(
        X = 1:job_details$n_jobs,
        FUN = map_SKAT_nm,
        this_phenotype = this_phenotype,
        covariates = covariates,
        pos_and_SNP_list = new_pre_allocated_SNP_windows,
        window_size = window_size,
        raw_file_path = new_pre_allocated_SNP_windows[[1]][[3]],
        n_permutations = job_details$n_permutations_per_job,
        resampling=TRUE)

      p_null_tallies <- dplyr::bind_rows(p_null_tallies)

      message(paste0(Sys.time(), " - Finished resampling up to ",
                     job_details$actual_n_permutations, " permutations in",
                     (proc.time() - timer)[3], "s\n\n"))

      p_empirical_table <- p_empirical_from_tally(
        p_null_tallies = p_null_tallies,
        scaffold_ID = basename(raw_file_path))

    }

    master_output <- rbind(
      master_output[-which(
        master_output$position %in% p_empirical_table$position), ],
      p_empirical_table)

  }

  output_dir <- paste0(output_dir, "/", job_id)
  output_basename <- paste0(output_dir, "SKAT_finaltraining_Chr",
                            basename(raw_file_path), "_Pheno",
                            basename(phenodata), ".csv")

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
