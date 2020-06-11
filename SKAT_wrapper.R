# SKAT_mydata.R

#install.packages("SKAT")
#install.packages("tari")
library(SKAT)
library(data.table)
library(tools)
library(foreach)
library(doParallel)

library(optparse)


# Read arguments from command line  ---------------------------------------


option_list = list(
  make_option(c("-p", "--phenodata_path"),
              type="character",
              default=NULL,
              help="path to phenotype file",
              metavar="character"),
  make_option(c("-c", "--covariates"),
              type="character",
              default=NULL,
              help="path to covariates file",
              metavar="character"),
  make_option(c("-r", "--raw_file_path"),
              type="character",
              default=NULL,
              help="path to genodata in raw format (see PLINK manual)",
              metavar="character"),
  make_option(c("-S", "--window_size"),
              type="numeric",
              default=3000,
              #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
              metavar="numeric"),
  make_option(c("-s", "--window_shift"),
              type="numeric",
              default=1000,
              #help="Number of pixels passing intensity threshold for binary classifications of transgenic or not",
              metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


setwd(output_dir)
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

# Build null models up to 1M ----------------------------------------------

#null_model_noresample <- SKAT_Null_Model(this_phenotype ~ 1 + fam_cov$V3 + fam_cov$V4 + fam_cov$V5 + covariates$Phase1 + covariates$Phase2 + covariates$Phase3 + covariates$Diameter, out_type="C")  
null_model_1kresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                         n.Resampling=1000, type.Resampling="bootstrap")  
null_model_10kresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                          n.Resampling=10000, type.Resampling="bootstrap")  
null_model_100kresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                           n.Resampling=100000, type.Resampling="bootstrap")  
null_model_1Mresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                         n.Resampling=1000000, type.Resampling="bootstrap")  

runSKAToneChr <- function(raw_file_path, window_size, window_shift, output_dir = "/scratch2/NSF_GWAS/Results/SKAT/"){

    master_output <- data.frame(matrix(NA, ncol=6))
    colnames(master_output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled", "SKAT_O_p-val", "SKAT_O_p-val_resampled")
    master_output$`SKAT_p-val_resampled` <- as.numeric(master_output$`SKAT_p-val_resampled`)
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
    
    #print(paste0("This position is ", this_position))
    #print(paste0("Max on this scaffold is ", max(this_scaff_subset$POS)))
    maxwindow <- max(this_scaff_subset$POS)
    while(this_position < maxwindow){
        #browser()

        print(paste0("Running window at position ",
                     this_position,
                     " out of ", maxwindow))
        #gc() #slowdown massive
        this_window_start <- proc.time()
      
        elapsed_time <- ((proc.time() - whole_genome_start)[3])/3600
        elapsed_time_to_print <- prettyNum(elapsed_time, digits = 2)
        so_far_n_test <- so_far_n_test + 1
        percent_tests_complete <- so_far_n_test/total_n_test
        time_per_test <- elapsed_time/so_far_n_test
        to_do_tests <- total_n_test-so_far_n_test
        estimated_time_remaining <- (to_do_tests*time_per_test)
        estimated_time_remaining_to_print <- prettyNum(estimated_time_remaining, digits = 2)
        #print(paste0("Running for: ", elapsed_time_to_print, "hr"))
        print(paste0("Remaining time: ", estimated_time_remaining_to_print, "hr"))
        print(paste0(prettyNum((percent_tests_complete*100), digits = 2),
                     "% complete for phenotype ",
                     basename(file_path_sans_ext(file_path_sans_ext(opt$phenodata_path)))))        
      
      
        locus_of_interest <- this_position
        window_start <- locus_of_interest - (window_size/2)
        
        if (window_start < 0){
          window_start <- 0
        }
        window_end <- locus_of_interest + (window_size/2)

        print(paste0("SNP window: ", window_start, "-", window_end))

        indices_to_pull <- which(with(this_scaff_subset, this_scaff_subset$POS >= window_start & this_scaff_subset$POS <= window_end), arr.ind = TRUE)
        if (length(indices_to_pull) == 0){
          this_position <- this_position + window_shift
          next
        }
        genodata_thiswindow <- this_scaff_subset[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)],]
        genodata_thiswindow[,1:6] <- NULL
        genodata_thiswindow <- data.frame(genodata_thiswindow)
        Z <- t(as.matrix(genodata_thiswindow))
        colnames(Z) <- NULL
        
        elapsed_time <- prettyNum(((proc.time() - this_window_start)[3])/60, digits = 2)
        
        if((1+1) == 2){
          this_SKAT_out <- SKAT(Z, null_model_1kresample, kernel = "linear.weighted")
          p_resampled_SKAT <- (length(subset(this_SKAT_out$p.value.resampling,
                                             this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))/1000
          
          if(this_SKAT_out$p.value == 1){
            this_position <- this_position + window_shift
            list_positions_pval_1 <- c(list_positions_pval_1, this_position)
            cat("\n")
            next
          }
          if(p_resampled_SKAT <= 0.010){
            print(paste0("P-value ", p_resampled_SKAT, " below 0.01 with resampling, so re-do test with 10,000 resample null model"))
            this_SKAT_out <- SKAT(Z, null_model_10kresample, kernel = "linear.weighted")
            p_resampled_SKAT <- (length(subset(this_SKAT_out$p.value.resampling,
                                               this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))/10000
            ##print(paste0("Head of resampled p-values: ", unlist(head(this_SKAT_out$p.value.resampling))))
            if(p_resampled_SKAT <= 0.0010){
              print(paste0("P-value ", p_resampled_SKAT, " below 0.001 with resampling, so re-do test with 100,000 resample null model"))
              this_SKAT_out <- SKAT(Z, null_model_100kresample, kernel = "linear.weighted")
              p_resampled_SKAT <- (length(subset(this_SKAT_out$p.value.resampling,
                                                 this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))/100000
              ##print(paste0("Head of resampled p-values: ", unlist(head(this_SKAT_out$p.value.resampling))))
              if(p_resampled_SKAT <= 0.00010){
                print(paste0("P-value ", p_resampled_SKAT, " below 0.0001 with resampling, so re-do test with 1,000,000 resample null model"))
                this_SKAT_out <- SKAT(Z, null_model_1Mresample, kernel = "linear.weighted")
                p_resampled_SKAT <- (length(subset(this_SKAT_out$p.value.resampling,
                                                   this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))/1000000
                ##print(paste0("Head of resampled p-values: ", unlist(head(this_SKAT_out$p.value.resampling))))
                if(p_resampled_SKAT <= 0.000010){
                  print(paste0("P-value ", p_resampled_SKAT, " below 0.00001 with resampling, so build 10 null models of 1M and re-do test"))
                  p.resample.all <- c()
                  start_10_nullmodel <- proc.time()
                  for(i in 1:10){
                    print(paste0("Null model ", basename(raw_file_path), " of 10"))
                    null_model_1Mresample <- NULL
                    gc()
                    set.seed(i*10)
                    null_model_1Mresample <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                                             n.Resampling=1000000, type.Resampling="bootstrap")
                    this_SKAT_out <- SKAT(Z, null_model_1Mresample, kernel = "linear.weighted")
                    print(paste0("Head of resampled p-values: ", unlist(head(this_SKAT_out$p.value.resampling))))
                    p.resample.all <- c(p.resample.all, this_SKAT_out$p.value.resampling)
                    print(paste0("Length of ALL resampled p-values for this kernel: ", length(p.resample.all)))
                  }
                  elapsed_time_resample10 <- ((proc.time() - start_10_nullmodel)[3])/60
                  elapsed_time_to_print <- prettyNum(elapsed_time_resample10, digits = 2)
                  n1<-length(which(this_SKAT_out$p.value >= p.resample.all)) +1 
                  n2<-length(p.resample.all)+1
                  p_resampled_SKAT <- n1/n2
                  print(paste0("Finish 10M test in ", elapsed_time_to_print, " min"))
                }
              }
            }
          }
        }
        print(paste0("This p-value, after resampling: ", p_resampled_SKAT))

        if(this_SKAT_out$p.value == 1){
           this_position <- this_position + window_shift
           list_positions_pval_1 <- c(list_positions_pval_1, this_position)
           cat("\n")
           next
         }
        print(paste0("This p-value, after resampling: ", p_resampled_SKAT))
        if (SKAT_O == "OFF"){
          to_append <- as.numeric(c(basename(raw_file_path), this_position, as.numeric(as.character(this_SKAT_out$p.value)), as.numeric(as.character(p_resampled_SKAT)), NA, NA))

        } else {
          to_append <- as.numeric(c(basename(raw_file_path), this_position, this_SKAT_out$p.value, p_resampled_SKAT, this_SKAT_O_out$p.value, p_resampled_SKAT_O))
        }
        master_output <- rbind(master_output, to_append)
        print("Lowest p-value in database:")
        print(min(na.omit(as.numeric(master_output$`SKAT_p-val_resampled`))))
        cat("\n")
        this_position <- this_position + window_shift
      }
      outputname <- paste0(output_dir, "SKAT_finaltraining_Chr", basename(raw_file_path), "_Pheno", basename(opt$phenodata_path), ".csv")
      print(paste0("Writing ", outputname))
      fwrite(master_output, outputname)
    
      proc.time() - whole_genome_start
}

# Run tests ---------------------------------------------------------------

#whole_genome_start <- proc.time()

runSKAToneChr(raw_file_path = opt$raw_file_path,
              window_size = opt$window_size,
              window_shift = opt$window_shift)



# Stop the clock
#proc.time() - whole_genome_start
