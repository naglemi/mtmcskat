MTSKAT_workflow <- function(phenodata,
                            covariates,
                            raw_file_path,
                            window_size,
                            window_shift,
                            output_dir,
                            pre_allocated_dir,
                            job_id,
                            RAM="AllRAM",
                            n_thread="AllCores"){ # get rid of switchpoint parameter once making function determine_switch_point?

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

  output_dir <- paste0(output_dir, "/", job_id)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  output_basename <- paste0(output_dir, "/pheno-",
                            basename(phenodata),
                            "-scaff-",
                            basename(raw_file_path))

  raw_out_name <- paste0(output_basename, ".csv")

  message(paste0("Writing ", raw_out_name))

  data.table::fwrite(master_output, raw_out_name)

  message(paste0("Done running in...",
                 print(proc.time() - whole_genome_start)[3],
                 "s"))

}
