# This is part of my latest attempt to parallelize at a high level over
#   scaffolds and phenotypes on SLURMS. We first made a callable script
#   (deploy_workflow.R), which will be called by THIS SCRIPT for each
#   phenotype/scaffold combination in a separate SLURMS batch job script.
#   This approach is described here, under "R Submission"
#   https://vsoch.github.io/lessons/sherlock-jobs/

submit_batch <- function(phenodata_list,
                         covariate,
                         raw_file_path_list,
                         window_size,
                         window_shift,
                         output_dir,
                         pre_allocated_dir,
                         #job_id,
                         desired_sig_figs,
                         min_accuracy,
                         max_accuracy,
                         plot,
                         RAM,
                         n_thread,
                         time,
                         email = FALSE,
                         allocation = FALSE){

  # Allow user to submit either a folder or a list of files for both
  #   phenodata and raw_file_path
  if(dir.exists(phenodata_list[1])){
    phenodata_list <- list.files(phenodata_list)
  }

  if(dir.exists(raw_file_path_list[1])){
    raw_file_path_list <- list.files(raw_file_path_list)
  }

  for(phenodata in phenodata_list){
    for(raw_file_path in raw_file_path_list){
      job_id <- paste0(basename(basename(raw_file_path)),
                       "_",
                       basename(basename(phenodata)))

      # Start writing to this file
      sink(paste0(".job/", job_id, '.job'))
      # the basic job submission script is a bash script
      cat("#!/bin/bash\n")
      cat("#SBATCH --job-name=", job_id, ".job\n", sep="")
      cat("#SBATCH --output=.out/", job_id, ".out\n", sep="")
      cat("#SBATCH --error=.out/", job_id, ".err\n", sep="")
      cat("#SBATCH --time=", time, "\n", sep="")

      if(allocation != FALSE){
        cat("#SBATCH -A")
      }

      if(RAM != "AllRAM"){
        cat("#SBATCH --mem=", RAM/1024, "\n", sep="")
      }

      if(email != FALSE){
        cat("#SBATCH --mail-type=ALL\n")
        cat("#SBATCH --mail-user=", email, "\n", sep="")
      }

      cat("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK")

      cat("Rscript mtmcskat_wrapper.R \\n",
          paste0("--phenodata=", phenodata, " \\n"),
          paste0("--covariate=", covariate, " \\n"),
          paste0("--raw_file_path=", raw_file_path, " \\n"),
          paste0("--window_size=", window_size, " \\n"),
          paste0("--window_shift=", window_shift, " \\n"),
          paste0("--output_dir=", output_dir, " \\n"),
          paste0("--pre_allocated_dir=", pre_allocated_dir, " \\n"),
          paste0("--job_id=", job_id, " \\n"),
          paste0("--desired_sig_figs=", desired_sig_figs, " \\n"),
          paste0("--min_accuracy=", min_accuracy, " \\n"),
          paste0("--max_accuracy=", max_accuracy, " \\n"),
          paste0("--plot=", plot, " \\n"),
          paste0("--RAM=", RAM, " \\n"),
          paste0("--n_thread=", n_thread, " \\n"),
          sep="")

      # Close the sink!
      sink()

      # Submit to run on cluster
      system(paste("sbatch",
                   paste(".job/", job_id, ".job", sep="")))
    }
  }
}
