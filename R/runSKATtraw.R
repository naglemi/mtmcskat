runSKATtraw <- function(phenodata, covariates, raw_file_path, window_size, window_shift, output_dir = "/scratch2/NSF_GWAS/Results/SKAT/", ncore=24, SKAT_O = "OFF",
                          job_id, desired_sig_figs = 2, min_leading_0s = 2){

  # Prepare SKAT inputs other than null model -------------------------------

  test_info <- data.table(matrix("NA", ncol = 2, nrow = 1))
  colnames(test_info) <- c("Window size", "Elapsed time")
  list_positions_pval_1 <- c()
  SKAT_O <- "OFF"
  ptm <- proc.time()
  #this_position <- (window_size/2)

  this_scaff_subset <- fread(paste0(raw_file_path),
                             fill = TRUE)
  whole_genome_start <- proc.time()
  total_n_test <- ((round(max(this_scaff_subset$POS)/window_size))*3)

  # i <= levels(factor(this_scaff_subset$CHR))
  # if(length(i)>1){
  #   stop("There is more than one chromosome in this file. This code is meant for running chromosomes independently.")
  # }

  print(paste0("Read genotype data for raw file ", raw_file_path, " with dimensions ", dim(this_scaff_subset)))
  this_scaff_subset$POS <- as.numeric(as.character(this_scaff_subset$POS))

  n_test <- (round(max(this_scaff_subset$POS)/window_size))*3

  #print(paste0("Number of rows on this scaffold: ", nrow(this_scaff_subset)))

  #print(paste0("This position is ", this_position))
  #print(paste0("Max on this scaffold is ", max(this_scaff_subset$POS)))
  maxwindow <- max(this_scaff_subset$POS)
  #browser()

  # To accomodate for cases when the input file is a whole genome or choromosome files or just a part of a chromosome file,
  # the starting window position will depend on the minimum base position of any SNP in the SNP data.

  starting_position <- max((window_size/2),
                           round_any(min(this_scaff_subset$POS),
                                     window_size/2))

  window_list <- seq(starting_position,
                     maxwindow,
                     by=window_shift)


  grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
                      checkFunc = function(i) !is.na(i))

  #browser()
  cl <- makeCluster(ncore,
                    type = "FORK",
                    outfile="cluster.outfile")
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  #print("Begin NONparallel action")

  null_model_noresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C")

  master_output <- loopSKAT(grab_window = grab_window,
                            window_size = window_size,
                            raw_file_path = raw_file_path,
                            resampling = FALSE,
                            null_model = null_model_noresample,
                            n_permutations = NA)

  master_output <- post_process_master_output(master_output = master_output)

  #browser()
  max_leading_0s <- ceiling(-log(min(master_output$`SKAT_p-val`), base = 10)) - 1
  print(paste0("Max # leading 0s is: ", max_leading_0s))
  # Make sure this is right
  #browser()

  for(leading_0s in min_leading_0s:min(4, max_leading_0s)){
    print(paste0("Running for", leading_0s, " leading zeroes"))
    upper_bound_p_val_for_MC_subset <- 0.1^leading_0s
    lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
    print(paste0("Bounds from ", lower_bound_p_val_for_MC_subset, upper_bound_p_val_for_MC_subset))
    n_permutations <- 10^( leading_0s + desired_sig_figs )
    print(paste0("n_permutations: ", n_permutations))

    window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position
    if(length(window_list) == 0){
      print(paste0("No SNPs with p-vals within the range bounded from ", lower_bound_p_val_for_MC_subset, " to ", upper_bound_p_val_for_MC_subset))
      next
    }
    print("Making null model")
    null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                  n.Resampling=n_permutations, type.Resampling="bootstrap")
    print("Done making null model")
    grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
                        checkFunc = function(i) !is.na(i))
    if(leading_0s==5){
      browser()
    }
    output_resampled <- loopSKAT(grab_window = grab_window,
                                 window_size = window_size,
                                 raw_file_path = raw_file_path,
                                 resampling = TRUE,
                                 null_model = null_model,
                                 n_permutations = n_permutations)

    print("Dimensions of output from resampling:")
    print(dim(output_resampled))
    output_resampled <- post_process_master_output(master_output = output_resampled)
    #browser()
    master_output <- rbind(master_output,
                           output_resampled)
  }

  # NOW FOR THE 10M
  #browser()

  terminal_resampling <- TRUE

  large_n_tests <- n_permutations <- 10^(leading_0s+(desired_sig_figs))
  print(paste0("n_permutations: ", n_permutations))
  n_permutations_per_small_null_model = 100000
  n_small_null_models <- large_n_tests/n_permutations_per_small_null_model

  leading_0s <- 5

  if(terminal_resampling == TRUE){
    lower_bound_p_val_for_MC_subset <- 0
  }
  if(terminal_resampling == FALSE){
    lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
  }

  upper_bound_p_val_for_MC_subset <- 0.1^leading_0s



  window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position

  # Note that when parallelizing over null models, we define iterator inside loopSKAT.
  # unlike when parallelizing over windows where it is defined outside of the function.

  output_resampled <- loopSKAT(#grab_window = grab_window,
                               window_size = window_size,
                               this_scaff_subset = this_scaff_subset,
                               window_list = window_list,
                               raw_file_path = raw_file_path,
                               resampling = TRUE,
                               null_model = NA,
                               n_permutations = n_permutations_per_small_null_model,
                               n_small_null_models = n_small_null_models,
                               parallelize_over = "null_models")


  final_output <- rbind(master_output,
                        output_resampled)

  output_dir <- paste0(output_dir, "/", job_id)
  outputname <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), ".csv")
  print(paste0("Writing ", outputname))


  fwrite(final_output, outputname)

  output2name <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), "1M_pvals.csv")

  print("Done running in...")
  print(proc.time() - whole_genome_start)
}
