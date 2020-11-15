#' Run an automated MTMC-SKAT job
#'
#' This function provides a one-liner to run MTMC-SKAT over a scaffold.
#'
#' For whole-genome analysis, this function (accesible by command line) is to be parallelized over jobs submitted in a batch query system on a high-performance cluster, or looped over scaffolds if running on a single machine.
#'
#' @param phenodata string for filepath for phenotype file, with labeled columns for FID, IID, and trait.
#' @param covariates string for filepath for covariate file, with ordered covariate values in each column and no header
#' @param raw_file_path string for filepath to .traw file (see [PLINK file format reference](https://www.cog-genomics.org/plink/2.0/formats))
#' @param output_dir string for directory where results will be saved
#' @param n_thread numeric, maximum number of cores over which parallelization is allowed
#' @param job_id string, identifier to label the job; this identifier will go into output filename
#' @param desired_sig_figs Integer, minimum number of significant figures
#'   desired for empirical p-values
#' @param min_accuracy numeric, threshold x to begin resampling for SNP windows with p-values below 10^-x; default value is 2
#' @param max_accuracy numeric, limit x to end resampling for SNP windows with p-values below 10^-x, usually due to computational cost
#' @param switch_point numeric, limit x at which SNP windows with p-values below 10^-x must be testing by parallelizing over null models rather than SNP windows (may be set by user due to limitations on RAM preventing production of large null models)
#' @param plot if TRUE, produce Manhattan plot of results
#' @inheritParams pre_allocate
#' @inheritParams calculate_max_perm_per_core
#'
#' @return None; outputs are saved to the user-specified output directory, with the user-specified ```job_id```
#' @export
#'
#' @examples
#' \dontrun{
#' phenodata <- system.file("extdata",
#'   "TDZ_shoot_area.plink.pheno",
#'   package = "SKATMCMT")
#'
#' covariates <- system.file("extdata",
#'                           "poplar_PCs_covariates.tbt",
#'                           package = "SKATMCMT")
#'
#' raw_file_path <- system.file("extdata",
#'                              "poplar_SNPs_Chr10_14460to14550kb.traw",
#'                              package = "SKATMCMT")
#'
#' MTMCSKAT_workflow(phenodata = phenodata,
#'                   covariates = covariates,
#'                   raw_file_path = raw_file_path,
#'                   window_size = 3000,
#'                   window_shift = 1000,
#'                   output_dir = "Results/",
#'                   pre_allocated_dir = "pre_allocated_dir/",
#'                   n_thread = "AllCores",
#'                   job_id = "my_sample_analysis",
#'                   desired_sig_figs = 2,
#'                   min_accuracy = 2,
#'                   max_accuracy = 5,
#'                   plot = TRUE)
#' }
#'
MTMCSKAT_workflow <- function(phenodata, covariates, raw_file_path, window_size,
                              window_shift, output_dir, pre_allocated_dir,
                              job_id, desired_sig_figs = 2,
                              min_accuracy = 2,
                              max_accuracy = 5,
                              switch_point = 4, plot = TRUE,
                              RAM="AllRAM",
                              n_thread="AllCores"){ # get rid of switchpoint parameter once making function determine_switch_point?

  set_accuracy <- function(master_output,
                           max_accuracy){

    max_leading_0s <- ceiling(-log(
      min(stats::na.omit(master_output$`SKAT_p-val`)),
      base = 10)) - 1

    message(paste0("\nMax # leading 0s is: ", max_leading_0s, "\n"))

    if (max_accuracy != "Auto"){
      message(paste0("User-defined accuracy limit is 10^-", max_accuracy, "."))
      if (max_leading_0s < max_accuracy){
        message(paste("This level of accuracy may be more than the accuracy",
                      "needed to calculate empirical p-values since model",
                      "p-values from mtskat have no more than", max_leading_0s,
                      "leading zeros\n"))
      }
      if (max_leading_0s > max_accuracy){
        message(paste("This accuracy may not be enough to accurately",
                      "calculate all empirical p-values since p-values from",
                      "mtskat have as many as", max_leading_0s,
                      "leading zeros\n"))
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

  check_if_more_threads_than_windows <- function(
    pre_allocated_SNP_windows,
    n_thread){

    # If there are fewer SNP windows than cores, we should parallelize over nms
    if(length(pre_allocated_SNP_windows) < n_thread){
      # SWITCH MODE TO PARALLELIZE OVER NULL MODELS INTEAD OF WINDOWS

      message(paste0(
        "Automatically switching modes for windows with p-vals <= 10^-",
        leading_0s,
        " because there are fewer windows (",
        length(pre_allocated_SNP_windows),
        ") than there are cores in use (",
        n_thread,
        ").\n"))
      return(TRUE)
    }

    if(length(pre_allocated_SNP_windows) >= n_thread){
      return(FALSE)
    }

  }

  check_if_null_models_fit_in_RAM_per_thread <- function(n_permutations,
                                                       RAM_per_permutation,
                                                       RAM_per_thread){

    # to parallelize over SNP windows, we need to be able to fit a copy
    # of the (same) null model for each thread
    RAM_per_thread_needed_for_mtmcskat_SNPs <-
      RAM_per_permutation * n_permutations

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

  this_phenotype <- unlist(utils::read.csv(phenodata,
                                           sep = "\t",
                                           colClasses = c("character",
                                                          "character",
                                                          "numeric"))[, 3])

  # fread is robust, automatically detects if there is a header and type of
  # separation, which may vary across covariate files (lack a standard format).
  covariates <- data.table::fread(covariates)

  if(n_thread=="AllCores") {
    n_thread <- future::availableCores()
  }
  if(RAM=="AllRAM") {
    RAM <- as.numeric(benchmarkme::get_ram())
  }

  whole_genome_start <- proc.time()

  future::plan("multisession")

  # Initial SKAT round w/o resampling, parallelized over SNP windows --------

  pre_allocated_SNP_windows <-
    pre_allocate(
      raw_file_path = raw_file_path,
      window_size = window_size,
      window_shift = window_shift,
      pre_allocated_dir = pre_allocated_dir)

  master_output <-
    mtskat(
      this_phenotype = this_phenotype,
      covariates = covariates,
      raw_file_path = pre_allocated_SNP_windows[[1]][[3]],
      pre_allocated_SNP_windows = pre_allocated_SNP_windows,
      n_thread = n_thread)

  # Secondary SKAT round w/ resampling, multithread over SNP windows or NMs-----

  max_accuracy <-
    set_accuracy(
      master_output = master_output,
      max_accuracy = max_accuracy)

  ## A null model with HOW MUCH resampling can fit into the RAM allotted
  ## for a single thread?

  RAM_per_permutation <-
    benchmark_nm(
      phenotype = this_phenotype,
      covariates = covariates,
      benchmark_size = 1000)

  # This parameter only used if multithreading over null models
  max_permutations_per_job <-
    calculate_max_perm_per_core(
      nm_RAM_per_perm = RAM_per_permutation,
      RAM = size_RAM_wiggle(RAM = RAM, wiggle_factor = 4),
      n_thread = n_thread)

  RAM_per_thread <- RAM / n_thread

  for(leading_0s in min_accuracy:max_accuracy){

    terminal_resampling <-
      determine_whether_final_resampling(
        leading_0s = leading_0s,
        max_accuracy = max_accuracy)

    boundaries <-
      determine_subset_pval_boundaries(
        leading_0s = leading_0s,
        desired_sig_figs = desired_sig_figs,
        terminal_resampling = terminal_resampling)

    new_pre_allocated_SNP_windows <-
      re_allocate_windows(
        x = master_output,
        upper_bound = boundaries$upper,
        lower_bound = boundaries$lower,
        pre_allocated_SNP_windows = pre_allocated_SNP_windows)

    if(new_pre_allocated_SNP_windows == "None") next

    n_permutations <-
      determine_n_permutations(
      leading_0s = leading_0s,
      desired_sig_figs = desired_sig_figs)

    null_models_fit_in_RAM_per_thread <-
      check_if_null_models_fit_in_RAM_per_thread(
        n_permutations = n_permutations,
        RAM_per_permutation = RAM_per_permutation,
        RAM_per_thread = RAM_per_thread)

    more_threads_than_windows <-
      check_if_more_threads_than_windows(
        pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
        n_thread = n_thread)

    message(paste0(Sys.time(), " - Running resampling for ",
                   length(new_pre_allocated_SNP_windows),
                   " SNP windows with p-values between ",
                   boundaries$upper, " and ",
                   boundaries$lower,
                   " (", leading_0s, " leading zeroes) with at least ",
                   n_permutations, " permutations so empirical p-values",
                   " can be calculated to ", desired_sig_figs,
                   " significant figures...\n"))

    if(null_models_fit_in_RAM_per_thread & !more_threads_than_windows){
      message("Multithreading over SNP windows\n")

      add_to_master_output <-
        mtmcskat_SNPs(
          this_phenotype = this_phenotype,
          covariates = covariates,
          n_permutations = n_permutations,
          pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
          scaffold_ID = pre_allocated_SNP_windows[[1]][[3]],
          n_thread = n_thread)
    }

    if(!null_models_fit_in_RAM_per_thread | more_threads_than_windows){
      message("Multithreading over null models\n")

      add_to_master_output <-
        mtmcskat_NullModels(
          this_phenotype = this_phenotype,
          covariates = covariates,
          n_thread = n_thread,
          n_permutations = n_permutations,
          max_permutations_per_job = max_permutations_per_job,
          pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
          scaffold_ID = pre_allocated_SNP_windows[[1]][[3]])

    }

    # Update the master output, replacing entries for which we now have...
    # ...empirical p-values
    master_output <-
      rbind(
        master_output[-which(
          master_output$position %in% add_to_master_output$position), ],
        add_to_master_output)

  }

  output_dir <- paste0(output_dir, "/", job_id)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_basename <- paste0(output_dir, "/pheno-",
                            basename(phenodata),
                            "-scaff-",
                            basename(raw_file_path))

  raw_out_name <- paste0(output_basename, ".csv")
  plot_out_name <- paste0(output_basename, ".png")

  message(paste0("Writing ", raw_out_name))

  data.table::fwrite(master_output, raw_out_name)

  if(plot==TRUE){

    compare_plots(master_output,
                  plot_out_name = plot_out_name,
                  scaffold_ID = pre_allocated_SNP_windows[[1]][[3]])

  }

  message(paste0("Done running in...",
                 print(proc.time() - whole_genome_start)[3],
                 "s"))

}
