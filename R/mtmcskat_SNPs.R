mtmcskat_SNPs <- function(pre_allocated_SNP_windows,
                          n_permutations,
                          this_phenotype,
                          covariates,
                          scaffold_ID){

  null_model <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates),
    n.Resampling = n_permutations,
    type.Resampling = "bootstrap")

  time_to_run_mapping <- proc.time()

  add_to_master_output <- future.apply::future_lapply(
    X = pre_allocated_SNP_windows,
    FUN = mappable_SKAT,
    scaffold_ID = scaffold_ID,
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

        n_permutations_per_job <- max_permutations_per_job
        n_jobs <- ceiling(n_permutations / n_permutations_per_job)
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
                 job_details$n_permutations_per_job, " permutations ",
                 "for a total of ", job_details$actual_n_permutations,
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

  message(paste0(Sys.time(), " - Finished resampling up to ",
                 job_details$actual_n_permutations, " permutations in",
                 (proc.time() - timer)[3], "s\n\n"))

  p_empirical_table <- p_empirical_from_tally(
    p_null_tallies = p_null_tallies,
    scaffold_ID = scaffold_ID)

  p_empirical_table[order(p_empirical_table$position), ]

}
