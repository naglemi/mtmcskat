runSKAToneChr <- function(phenodata, covariates, raw_file_path, window_size, window_shift, output_dir = "/scratch2/NSF_GWAS/Results/SKAT/", ncore=24, SKAT_O = "OFF",
                          job_id, desired_sig_figs = 2){

  # Prepare SKAT inputs other than null model -------------------------------

  test_info <- data.table(matrix("NA", ncol = 2, nrow = 1))
  colnames(test_info) <- c("Window size", "Elapsed time")
  list_positions_pval_1 <- c()
  SKAT_O <- "OFF"
  ptm <- proc.time()
  this_position <- (window_size/2)

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

  print(paste0("This position is ", this_position))
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

  # Make an iterator function so cores pull windows without copying  --------
  ## without copying whole traw file
  #grab_window <- iter(function() this_scaff_subset[which(with(this_scaff_subset, this_scaff_subset$POS >= (this_position - (window_size/2)) & (this_position +(window_size/2)), arr.ind = TRUE)])

  # grab_window <- iter(window_list,
  #                     checkFunc = function(i) extract_window(i, window_size = window_size, this_scaff_subset = this_scaff_subset)) # https://datawookie.netlify.app/blog/2013/11/iterators-in-r/

  # grab_window <- iter(function(window_list) extract_window(window_list, window_size = window_size, this_scaff_subset = this_scaff_subset))

  # grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
  #                     checkFunc = function(i) !is.na(i))

  grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
                      checkFunc = function(i) !is.na(i))

  #browser()
  # cl <- makeCluster(ncore,
  #                   type = "FORK",
  #                   outfile="cluster.outfile")
  # registerDoParallel(cl)
  # on.exit(stopCluster(cl))

  #print("Begin NONparallel action")

  null_model_noresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C")

  master_output <- loopSKAT(grab_window = grab_window,
                            window_size = window_size,
                            raw_file_path = raw_file_path,
                            resampling = FALSE,
                            null_model = null_model_noresample)

  #browser()
  # print("Begin NONparallel action")
  # master_output <- foreach(i=window_list, .combine="rbind", .export=c("window_size")) %do% {
  #   print("Getting Z frmo nextElem")
  #   Z <- nextElem(grab_window)
  #   # FOR SOME REASON WE BREAK OUT OF THE LOOP HERE!!! FIX THIS
  #   print(Z)
  #   stop()
  #   print("Got Z")
  #
  #   # The above function returns None if there are no SNPs within the window. Because of this,
  #   # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
  #   if(is.matrix(Z)==TRUE){
  #     print("About to go into SKAT_one_window")
  #     this_SKAT_out <- SKAT_one_window(this_position = i,
  #                                      window_size,
  #                                      Z = Z,
  #                                      raw_file_path = raw_file_path,
  #                                      resampling = FALSE)
  #     print("Out of SKAT_one_window")
  #   } else {
  #     print("Not matrix!")
  #   }
  #   print("Dimensions of SKAT out:")
  #   return(this_SKAT_out)
  # }
  post_process_master_output <- function(master_output){
    master_output <- data.frame(master_output)
    colnames(master_output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled", "SKAT_O_p-val", "SKAT_O_p-val_resampled")
    master_output$`SKAT_p-val_resampled` <- as.numeric(as.character(master_output$`SKAT_p-val_resampled`))
    master_output
  }

  master_output <- post_process_master_output(master_output = master_output)

  max_leading_0s <- attr(regexpr("(?<=\\.)0+",
                                 min(master_output$`SKAT_p-val`),
                                 perl = TRUE), "match.length") # https://stackoverflow.com/questions/35553244/count-leading-zeros-between-the-decimal-point-and-first-nonzero-digit
  print(paste0("Max # leading 0s is: ", max_leading_0s))

  for(leading_0s in 2:min(5, max_leading_0s)){
    print(paste0("Running for", leading_0s, " leading zeroes"))
    upper_bound_p_val_for_MC_subset <- 0.1^leading_0s
    lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
    n_permutations <- 10^(leading_0s+(desired_sig_figs-1))

    window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position
    print("Making null model")
    null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                  n.Resampling=n_permutations, type.Resampling="bootstrap")
    print("Done making null model")
    grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
                        checkFunc = function(i) !is.na(i))

    output_resampled <- loopSKAT(grab_window = grab_window,
                                 window_size = window_size,
                                 raw_file_path = raw_file_path,
                                 resampling = TRUE,
                                 null_model = null_model,
                                 n_permutations = n_permutations)

    print("Dimensions of output from resampling:")
    print(dim(output_resampled))
    output_resampled <- post_process_master_output(master_output = output_resampled)
    master_output <- rbind(master_output, output_resampled)
  }

  # NOW FOR THE 10M
  browser()
  grab_window <- iter(window_list_10M_resample,
                      checkFunc = function(i) extract_window(i, window_size = window_size, this_scaff_subset = this_scaff_subset)) # https://datawookie.netlify.app/blog/2013/11/iterators-in-r/

  output_10M_resample <- foreach(i=1:10, .combine="cbind") %do% {
    set.seed(i*10)
    null_model_1Mresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                             n.Resampling=1000000, type.Resampling="bootstrap")

    foreach(i=window_list_1M_resample, .combine="rbind", .export=c("window_size")) %dopar% {
      Z <- nextElem(grab_window)
      # The above function returns None if there are no SNPs within the window. Because of this,
      # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.

      # I need to massively rework this so we only generate 10 1M resample null models in total, and use the same ones for all windows.
      # This is a very special case because we have no machine with ~1TB RAM that would be needed for a 10M resample null model

      if(is.matrix(Z)==TRUE){
        SKAT_one_window(this_position = window,
                        window_size,
                        Z = Z,
                        raw_file_path = raw_file_path,
                        resampling = TRUE,
                        null_model = null_model_1Mresample,
                        n_permutations = 1000000,
                        return_all_p_vals = TRUE)
      }

      # THE BELOW LINE IS WRONG. THIS SHOULD BE A LINE THAT ADDS EACH SET OF 1M p-vals TO THE TOTAL
      # p.resample.all <- c(p.resample.all, this_SKAT_out$p.value.resampling)

      ## SECOND THOUGHTS:
      ## THINK I SHOULD PARALLELIZE OVER WINDOWS STILL, AND USE A SIMPLE FOR-LOOP for EACH 1M NULL MODEL GENERATED
      ## SINCE THEY ONLY TAKE A MINUTE OR SO EACH
      ## ... AND IT's EASIER TO CODE

    }

  }


  output_no_resample_subset <- master_output[which(master_output$SKAT_p-val > 0.01),]

  final_output <- rbind(output_no_resample_subset,
                        output_1k_resample,
                        output_10k_resample,
                        output_100k_resample,
                        output_10M_resample)

  output_dir <- paste0(output_dir, "/", job_id)
  outputname <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), ".csv")
  print(paste0("Writing ", outputname))


  fwrite(final_output, outputname)

  output2name <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), "1M_pvals.csv")

  proc.time() - whole_genome_start
}
