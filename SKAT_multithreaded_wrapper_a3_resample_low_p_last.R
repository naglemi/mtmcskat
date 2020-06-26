# SKAT_mydata.R
rm(list=ls())

#install.packages("SKAT")
#install.packages("tari")
library(SKAT)
library(data.table)
library(tools)
library(foreach)
library(doParallel)

library(optparse)
library(iterators)

# Read arguments from command line  ---------------------------------------


# option_list = list(
#   make_option(c("-p", "--phenodata_path"),
#               type="character",
#               default=NULL,
#               help="path to phenotype file",
#               metavar="character"),
#   make_option(c("-c", "--covariates"),
#               type="character",
#               default=NULL,
#               help="path to covariates file",
#               metavar="character"),
#   make_option(c("-r", "--raw_file_path"),
#               type="character",
#               default=NULL,
#               help="path to genodata in raw format (see PLINK manual)",
#               metavar="character"),
#   make_option(c("-S", "--window_size"),
#               type="numeric",
#               default=3000,
#               #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
#               metavar="numeric"),
#   make_option(c("-s", "--window_shift"),
#               type="numeric",
#               default=1000,
#               #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
#               metavar="numeric"),
#   make_option(c("-o", "--output_dir"),
#               type="character",
#               default="/scratch2/NSF_GWAS/Results/SKAT/",
#               help="Path to directory where results will be saved",
#               metavar="character"),
#   make_option(c("-j", "--job_id"),
#               type="character",
#               default="no_job_id",
#               help="User-provided job identifier, to be used as subfolder in which results are placed",
#               metavar="character")
# );
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

opt <- list(phenodata_path = "/scratch2/NSF_GWAS/phenodata/third_training/PCs_ratios_binarized_4class/TDZ_shoot_area.emmax.pheno",
            #phenodata_path = "/scratch2/NSF_GWAS/phenodata/epicormic/EpicormicBudBreak_South_March_2014.header.pheno",
            covariates = "/scratch2/NSF_GWAS/phenodata/Covariates_DiamFinal_PhaseFinal_PCs.txt",
            #raw_file_path = "/scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr10.traw",
            raw_file_path = "/scratch2/NSF_GWAS/genodata/MAF0.0_geno882_Chr10_portion.traw",
            window_size = 3000,
            window_shift = 1000,
            output_dir = "/scratch2/NSF_GWAS/Results/SKAT/",
            job_id = "fixing_iter_a3")

setwd(opt$output_dir)
#phenodata_path <- "/scratch2/NSF_GWAS/phenodata/final_training/shoot_5w.header.pheno"
covariates <- fread(opt$covariates)

# Adding new phase data 3/29/20 ------------------------------------------------
phenodata <- fread(opt$phenodata_path)
this_phenotype <- as.numeric(unlist(phenodata[,3]))

# Count SNP windows -------------------------------------------------------

#total_n_test <- 390000

so_far_n_test <- 0

# for (that_scaff in levels(factor(trawA_in$CHR))[1:19]){
#   #print(that_scaff)
#   that_scaff_subset <- subset(trawA_in, trawA_in$CHR == that_scaff)
#   if(nrow(that_scaff_subset) != 0){
#     total_n_test <- total_n_test + ((round(max(that_scaff_subset$POS)/window_size))*3)
#   }
#   #print(total_n_test)
# }



extract_window <- function(this_position, window_size, this_scaff_subset){
  print("Just entered extract_window")
  print(paste0("This position (within extract_window) is: ", this_position, " with SNPs: "))
  
  if(is.na(this_position)){
    print("NA position. Reached end of scaffold?")
    stop("StopIteration")
  }
  
  locus_of_interest <- this_position
  window_start <- locus_of_interest - (window_size/2)
  
  if (window_start < 0){
    window_start <- 0
  }
  window_end <- locus_of_interest + (window_size/2)
  
  print(paste0("SNP window: ", window_start, "-", window_end))
  
  indices_to_pull <- which(with(this_scaff_subset, this_scaff_subset$POS >= window_start & this_scaff_subset$POS <= window_end), arr.ind = TRUE)
  print("Got indices if they exist")
  if (length(indices_to_pull) == 0){
    #this_position <- this_position + window_shift
    #next
    #return()
    # I think we need return(NA) instead of return() to stop is.matrix(Z) from giving a "missing value where TRUE/FALSE needed" error
    return(NA)
  }
  print("Determined whether there are any indices")
  genodata_thiswindow <- this_scaff_subset[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)],]
  #return(genodata_thiswindow)
  print("Obtained subseted genodata")
  
  genodata_thiswindow[,1:6] <- NULL
  genodata_thiswindow <- data.frame(genodata_thiswindow)
  print("Reformatted genodata (pt1)")

  Z <- t(as.matrix(genodata_thiswindow))
  colnames(Z) <- NULL
  print("Reformatted genodata (pt2)")
  
  print("Z extracted dimensions and first 10 col, row:")
  print(dim(Z))
  print(head(Z)[,1:min(10, ncol(Z))])
  print("Now about to return Z")
  return(list(this_position, Z))
}

calculate_SKAT_empirical_p_old <- function(n_permutations, null_model, return_all_p){
  if(n_permutations < 1000000){
    this_SKAT_out <- SKAT(Z, null_model, kernel = "linear.weighted")
    p_resampled_SKAT <- ((length(subset(this_SKAT_out$p.value.resampling,
                                        this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))+1)/(n_permutations+1)
  } else if(n_permutations == 1000000){
    # print(paste0("P-value ", p_resampled_SKAT, " below 0.00001 with resampling, so build 10 null models of 1M and re-do test"))
    p.resample.all <- c()
    start_10_nullmodel <- proc.time()
    for(i in 1:10){
      # print(paste0("Null model ", basename(raw_file_path), " of 10"))
      null_model_1Mresample <- NULL
      gc()
      set.seed(i*10)
      null_model_1Mresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                               n.Resampling=1000000, type.Resampling="bootstrap")
      this_SKAT_out <- SKAT(Z, null_model_1Mresample, kernel = "linear.weighted")
      # print(paste0("Head of resampled p-values: ", unlist(head(this_SKAT_out$p.value.resampling))))
      p.resample.all <- c(p.resample.all, this_SKAT_out$p.value.resampling)
      # print(paste0("Length of ALL resampled p-values for this kernel: ", length(p.resample.all)))
    }
    # elapsed_time_resample10 <- ((proc.time() - start_10_nullmodel)[3])/60
    # elapsed_time_to_print <- prettyNum(elapsed_time_resample10, digits = 2)
    n1<-length(which(this_SKAT_out$p.value >= p.resample.all)) +1 
    n2<-length(p.resample.all)+1
    p_resampled_SKAT <- n1/n2
    # print(paste0("Finish 10M test in ", elapsed_time_to_print, " min"))
  } else {
    stop(paste0("Unsupported value for n_permutations ", n_permutations))
  }
  
}

calculate_SKAT_empirical_p <- function(n_permutations, null_model, return_all_p){
  this_SKAT_out <- SKAT(Z, null_model, kernel = "linear.weighted")
  if(return_all_p == FALSE){
    p_resampled_SKAT <- ((length(subset(this_SKAT_out$p.value.resampling,
                                        this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))+1)/(n_permutations+1)
    return(p_resampled_SKAT)
  }
  if(return_all_p == TRUE){
    return(this_SKAT_out$p.value.resampling)
  }
}

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
    window_list <- seq((window_size/2),
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
    
    loopSKAT <- function(grab_window, window_size, raw_file_path, resampling, null_model, ...){
      master_output <- foreach(pos_and_SNPs=grab_window, 
                               .combine="rbind",
                               .export=c("window_size"),
                               .errorhandling="pass") %do% {
                                 print("Now inside foreach loop. Z is: ")
                                 #print(Z)
                                 #stop()
                                 #print("Got Z")
                                 
                                 if(is.list(pos_and_SNPs) == FALSE){
                                   print("No SNPs in this window")
                                   return(NA)
                                 }
                                 
                                 print(pos_and_SNPs[[1]])
                                 
                                 # The above function returns None if there are no SNPs within the window. Because of this,
                                 # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
                                 if(is.matrix(pos_and_SNPs[[2]])==TRUE){
                                   print("About to go into SKAT_one_window")
                                   this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNPs[[1]],
                                                                    window_size,
                                                                    Z = pos_and_SNPs[[2]],
                                                                    raw_file_path = raw_file_path,
                                                                    resampling = resampling,
                                                                    null_model = null_model_noresample)
                                   print("Out of SKAT_one_window")
                                 } else {
                                   print("Not matrix!")
                                   return(NA)
                                 }
                                 print("length of SKAT out:")
                                 print(length(this_SKAT_out))
                                 print("This SKAT out about to be returned:")
                                 print(this_SKAT_out)
                                 return(this_SKAT_out)
                               }
      master_output
    }
    
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
    
    for(leading_0s in 2:5){
      upper_bound_p_val_for_MC_subset <- 0.1^leading_0s
      lower_bound_p_val_for_MC_subset <- 0.1^(leading_0s+(desired_sig_figs-1))
      n_permutations <- 10^(leading_0s+(desired_sig_figs-1))
      
      window_list <- master_output[which(master_output$`SKAT_p-val` <= upper_bound_p_val_for_MC_subset & master_output$`SKAT_p-val` > lower_bound_p_val_for_MC_subset),]$position
      null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                    n.Resampling=n_permutations, type.Resampling="bootstrap")  
      grab_window <- iter(obj = function(i) extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset),
                          checkFunc = function(i) !is.na(i))
      
      output_resampled <- loopSKAT(grab_window = grab_window,
                                   window_size = window_size,
                                   raw_file_path = raw_file_path,
                                   resampling = TRUE,
                                   null_model = null_model,
                                   n_permutations = n_permutations)
      
      
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

SKAT_one_window <- function(this_position, window_size, Z, SKAT_O = "OFF", raw_file_path, null_model, n_permutations, resampling=FALSE, return_all_p_vals=FALSE){
  #print("Inside SKAT_one_window")
  #browser()
  this_SKAT_out <- SKAT(Z, null_model, kernel = "linear.weighted")
  #print("Done running SKAT")
  if(resampling==TRUE){
    if(return_all_p_vals == FALSE){
      p_empirical <- calculate_SKAT_empirical_p(n_permutations = n_permutations,
                                                null_model = null_model)
      to_append <- as.numeric(c(basename(raw_file_path), this_position, as.numeric(as.character(this_SKAT_out$p.value)), p_empirical, NA, NA))
    }
    if(return_all_p_vals == TRUE){
      p_list <- calculate_SKAT_empirical_p(n_permutations = n_permutations,
                                           null_model = null_model,
                                           return_all_p = TRUE)
      to_append <- as.numeric(this_position, p_list)
    }

  }
  if(resampling==FALSE){
    #print("Not resampling (yet)")
    p_empirical <- -9
    to_append <- as.numeric(c(basename(raw_file_path), this_position, as.numeric(as.character(this_SKAT_out$p.value)), p_empirical, NA, NA))
  }
  
  print(length(to_append))
  print(to_append)
  
  to_append
}

whole_genome_start <- proc.time()

runSKAToneChr(phenodata = phenodata,
              covariates = covariates,
              raw_file_path = opt$raw_file_path,
              window_size = opt$window_size,
              window_shift = opt$window_shift,
              output_dir = opt$output_dir,
              job_id = opt$job_id)



# Stop the clock
proc.time() - whole_genome_start

