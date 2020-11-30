#' Run multi-threaded Monte Carlo Sequence Kernel Association Test over a
#' pre-allcoated set of SNP windows
#'
#' These function provide means of running SKAT over a set of pre-allocated
#' SNP windows (see \code{\link{pre_allocate}}) while performing resampling
#' with a user-specified number of permutations. Multithreading is implemented
#' over SNP windows or null models, depending on the function called.
#'
#' @inheritParams mtskat
#' @inheritParams chunk_windows
#' @inheritParams mappable_SKAT
#' @param max_permutations_per_job the maximum number of permutations that may fit
#' into a `SKAT_NULL_Model` object generated for any given thread
#'
#' @return A dataframe with four columns, for 1) scaffold ID, 2) SNP window
#' position, 3) p-values from the model used in SKAT without resampling, and
#' 4) empirical p-values
#' @export
#'
#' @details The choice of the appropriate function (either `mtmcskat_SNPs` or
#' `mtmcskat_NullModels`) depends on whether the user desired to multithread
#' over SNP windows or null models, which in turn depends on several factors.
#' First, if multi-threading over SNP windows, the number of SNP windows should
#' not be less than the number of threads available, as efficient computation
#' requires each thread receive at least one SNP window. Second, available RAM
#' limits the amount of permutation that can be accomodated simultaneously in
#' null models.
#'
#'
#' Multithreading over null models is recommended for situations
#' in which the number of permutations that must be tested cannot fit into
#' available RAM divided by available threads. Multithreading over SNP windows
#' is only recommended for situations involving relatively few permutations and
#' large number of SNP windows, for example, in the standard MTMCSKAT workflow,
#' when calculating empirical p-values for large groups of SNP windows with
#' initial \code{\link{mtskat}} p-values that are not very low.
#'
#' @examples
#'
#' data("sample_phenotype")
#' data("sample_covariates")
#' data("sample_pre_allocated_SNP_windows")
#'
#' # Multithreading over SNP windows
#' mtmcskat_SNPs(
#' this_phenotype = sample_phenotype,
#' covariates = sample_covariates,
#' n_permutations = 500,
#' pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
#' scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]],
#' n_thread = 2)
#'
#' # Multithreading over null models, where all necessary permutations can
#' #  simultaneously fit into memory and computation can be completed in a
#' #  single "batch."
#' mtmcskat_NullModels(
#' this_phenotype = sample_phenotype,
#' covariates = sample_covariates,
#' n_permutations = 500,
#' n_thread = 2,
#' max_permutations_per_job = 251,
#' pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
#' scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]])
#'
#'
#' # Multithreading over null models, where the the number of permutations
#' #   is greater than that which can fit into memory (as indicated by the
#' #   user or upstream functions through the `max_permutations_per_job`
#' #   argument), thus requiring multiple sequential "batches" of computation.
#' mtmcskat_NullModels(
#' this_phenotype = sample_phenotype,
#' covariates = sample_covariates,
#' n_permutations = 500,
#' n_thread = 2,
#' max_permutations_per_job = 249,
#' pre_allocated_SNP_windows = sample_pre_allocated_SNP_windows[2:4],
#' scaffold_ID = sample_pre_allocated_SNP_windows[[1]][[3]])
#'
mtmcskat_NullModels <- function(n_thread,
                                n_permutations,
                                max_permutations_per_job,
                                this_phenotype,
                                covariates,
                                pre_allocated_SNP_windows,
                                scaffold_ID){

  arrange_jobs_NullModel_multithreading <-
    function(n_thread,
             n_permutations,
             max_permutations_per_job){
      # How many jobs do we want to break this down into?
      #   Makes sense to have one for each core to minimize communication
      n_jobs <- n_thread

      # The number of permutations per null model should be such that
      #   there is only one round of communication... unless we don't have
      #   enough RAM for that, in which case we can divide jobs to give
      #   each core twice as many (or more if needed)

      n_total_permutations_per_thread <- ceiling(n_permutations / n_thread)

      if(n_total_permutations_per_thread >= max_permutations_per_job ){
        message(paste0("Memory constraint recognized; ",
                       "Unable to multithread over all null models ",
                       "simultaneously."))

        # Old way, not efficient... Doesn't use all cores for last round
        #   of resampling if there are fewer tasks remaining than n threads.
        # n_permutations_per_job <- max_permutations_per_job
        # n_jobs <- ceiling(n_permutations / n_permutations_per_job)

        # New approach, makes sure the number of tasks is always a multiple
        #   of the number of threads available, and decides on a number of
        #   permutations per core that allows us to use ALL threads ALWAYS
        #   to reach the desired # of permutations with EXACTLY the same number
        #   of permutations per thread

        # The smallest # jobs we could possibly divide needed permutations
        #   over without crashing due to running over RAM
        min_n_jobs <- ceiling(n_permutations / max_permutations_per_job)

        # To make full use of resources and ensure we're always running the
        #   maximum possible number of threads...
        n_jobs <- plyr::round_any(min_n_jobs, n_thread, f = ceiling)
        n_permutations_per_job <- ceiling(n_permutations / n_jobs)
      }

      if(n_total_permutations_per_thread < max_permutations_per_job ){
        message(paste0("Memory is enough to multithread over all null ",
                       "models simultaneously."))

        n_permutations_per_job <- n_total_permutations_per_thread
        n_jobs <- n_thread
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

  job_details <- arrange_jobs_NullModel_multithreading(
    n_thread = n_thread,
    n_permutations = n_permutations,
    max_permutations_per_job = max_permutations_per_job
  )

  message(paste0("To run ", job_details$n_jobs, "jobs, each with ",
                 format(job_details$n_permutations_per_job,
                        big.mark=",",scientific=FALSE), " permutations ",
                 "for a total of ",
                 format(job_details$actual_n_permutations,
                        big.mark=",",scientific=FALSE),
                 " permutations\n"))

  timer <- proc.time()

  p_null_tallies <- future.apply::future_lapply(
    X = 1:job_details$n_jobs,
    FUN = map_SKAT_nm,
    this_phenotype = this_phenotype,
    covariates = covariates,
    pos_and_SNP_list = pre_allocated_SNP_windows,
    scaffold_ID = pre_allocated_SNP_windows[[1]][[3]],
    n_permutations = job_details$n_permutations_per_job,
    resampling=TRUE)

  p_null_tallies <- dplyr::bind_rows(p_null_tallies)

  message(paste0(Sys.time(), " - Finished resampling with ",
                 format(job_details$actual_n_permutations,
                        big.mark=",",scientific=FALSE), " permutations in",
                 (proc.time() - timer)[3], "s\n\n"))

  p_empirical_table <- p_empirical_from_tally(
    p_null_tallies = p_null_tallies,
    scaffold_ID = scaffold_ID,
    min_n_permutations = n_permutations)

  p_empirical_table[order(p_empirical_table$position), ]

}


#' @rdname mtmcskat_NullModels
#' @export
mtmcskat_SNPs <- function(pre_allocated_SNP_windows,
                          n_permutations,
                          this_phenotype,
                          covariates,
                          scaffold_ID,
                          n_thread){

  message(paste(Sys.time(), "- Making null model with",
                n_permutations, "permutations..."))

  null_model <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates),
    n.Resampling = n_permutations,
    type.Resampling = "bootstrap")

  message(paste(Sys.time(), "- Complete\n"))

  pre_allocated_SNP_windows <- chunk_windows(
    pre_allocated_SNP_windows = pre_allocated_SNP_windows,
    n_thread = n_thread)

  time_to_run_mapping <- proc.time()

  #browser()

  add_to_master_output <- future.apply::future_lapply(
    X = pre_allocated_SNP_windows,
    FUN = mappable_SKAT,
    scaffold_ID = scaffold_ID,
    null_model = null_model,
    resampling = TRUE,
    n_permutations = n_permutations,
    chunk = TRUE)

  message(paste("Finished parallel run in",
                (proc.time() - time_to_run_mapping)[3],
                "seconds\n"))

  add_to_master_output <- dplyr::bind_rows(add_to_master_output)
  add_to_master_output <- post_process_master_output(
    master_output = add_to_master_output)

  add_to_master_output
}

