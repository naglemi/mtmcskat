loopSKAT <- function(grab_window, window_size, raw_file_path, resampling, null_model, n_permutations, parallelize_over = "windows", n_small_null_models = NA,
                     this_scaff_subset, window_list){
  print("Will use foreach now")
  if(parallelize_over=="windows"){
    print("Will parallelize over windows")
    master_output <- foreach(pos_and_SNPs=grab_window,
                             .combine="rbind",
                             .export=c("window_size",
                                       "n_permutations"),
                             .errorhandling="pass") %dopar% {
                               #print("Inside foreach loop")
                               #print("Now inside foreach loop. Z is: ")
                               ##print(Z)
                               #stop()
                               ##print("Got Z")

                               if(is.list(pos_and_SNPs) == FALSE){
                                 #print("No SNPs in this window")
                                 return(NA)
                               }

                               print(pos_and_SNPs[[1]])
                               print(dim(pos_and_SNPs[[2]]))

                               # The above function returns None if there are no SNPs within the window. Because of this,
                               # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
                               if(is.matrix(pos_and_SNPs[[2]])==TRUE){
                                 #print("About to go into SKAT_one_window")
                                 this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNPs[[1]],
                                                                  window_size,
                                                                  Z = pos_and_SNPs[[2]],
                                                                  raw_file_path = raw_file_path,
                                                                  resampling = resampling,
                                                                  null_model = null_model,
                                                                  n_permutations = n_permutations)
                                 #print("Out of SKAT_one_window")
                               } else {
                                 #print("Not matrix!")
                                 return(NA)
                               }
                               #print("length of SKAT out:")
                               #print(length(this_SKAT_out))
                               #print("This SKAT out about to be returned:")
                               #print(this_SKAT_out)
                               return(this_SKAT_out)
                             }
  }

  if(parallelize_over=="null_models"){
    print("Parallelize over null models")
    master_output <- foreach(null_model_index = 1:n_small_null_models,
                             .combine="rbind",
                             .export=c("window_size",
                                       "n_permutations"),
                             .errorhandling="pass") %dopar% {
                               print("Inside foreach loop")
                               #print(paste0("Now inside foreach loop for parallelizing over null models. This # null model is: ", null_model_index))
                               ##print(Z)
                               #stop()
                               ##print("Got Z")

                               #browser()
                               set.seed(null_model_index)
                               null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                                             n.Resampling=n_permutations, type.Resampling="bootstrap")
                               #print("Made null model")
                               ##browser()
                               master_output <- data.table(matrix(NA, nrow=1, ncol=(n_permutations+2)))
                               #brow


                               #print("Defined master_output and grab_window on this (to-be) node")
                               #print(paste0("Length of window_list is: ", length(window_list)))
                               for(i in 1:length(window_list)){
                                 ##print(paste0("On ", i, "th iteration through window_list"))
                                 pos_and_SNPs <- extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset)

                                 if(is.list(pos_and_SNPs) == FALSE){
                                   #print("No SNPs in this window")
                                   return(NA)
                                 }

                                 #print(pos_and_SNPs[[1]])

                                 # The above function returns None if there are no SNPs within the window. Because of this,
                                 # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
                                 if(is.matrix(pos_and_SNPs[[2]])==TRUE){
                                   #print("About to go into SKAT_one_window")
                                   this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNPs[[1]],
                                                                    window_size,
                                                                    Z = pos_and_SNPs[[2]],
                                                                    raw_file_path = raw_file_path,
                                                                    resampling = resampling,
                                                                    null_model = null_model,
                                                                    n_permutations = n_permutations,
                                                                    return_all_p_vals = TRUE)
                                   #print("Out of SKAT_one_window")

                                 } else {
                                   #print("Not matrix!")
                                   return(NA)
                                 }
                                 #print("length of SKAT out:")
                                 #print(length(this_SKAT_out))
                                 ##print("This SKAT out about to be returned:")
                                 ##print(this_SKAT_out)
                                 ##browser()
                                 master_output <- rbind(master_output, t(this_SKAT_out))
                                 ##print(paste0("Finishing ", i, "th iteration through window_list"))
                               }

                               #browser()
                               #master_output
                               # Get rid of the NA row used during initialization of the data frame
                               master_output <- master_output[-1,]
                               master_output <- t(master_output)
                               master_output <- as.data.frame(master_output)
                               rownames(master_output)[1:2] <- c("Position", "p_no_resampling")

                               abbreviated_master_output <- as.data.frame(t(master_output[1:2,]))
                               abbreviated_master_output$n_perm_p_below <- rep(NA, nrow(abbreviated_master_output))
                               abbreviated_master_output$n_perm_p_above <- rep(NA, nrow(abbreviated_master_output))

                               for(k in 1:ncol(master_output)){
                                 #print(paste0("For SNP at position ", master_output[1,k]))
                                 n_perm_p_below <- length(which(master_output[ 3:nrow(master_output) , k] <= master_output[ 2 , k ]))
                                 n_perm_p_above <- length(which(master_output[ 3:nrow(master_output) , k] > master_output[ 2 , k ]))
                                 abbreviated_master_output$n_perm_p_below[k] <- n_perm_p_below
                                 abbreviated_master_output$n_perm_p_above[k] <- n_perm_p_above
                               }

                               return(abbreviated_master_output)
                             }
    #browser()
    # At this point, we have a data structure containing a row for results from each X permutations results
    # With multiple rows for each SNP. Want to combine rows by SNP and calculate final empirical p-values.

    #master_output$n_perm_p_above <- NULL

    total_perm_p_below <- aggregate(master_output$n_perm_p_below,
                                    by=list(master_output$Position, master_output$p_no_resampling), FUN=sum) # https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
    total_perm_p_below$x <- total_perm_p_below$x + 1 # Add one for our p-value without resampling
    total_perm_p_below$empirical_p <- total_perm_p_below$x / ( n_permutations * n_small_null_models )

    master_output <- cbind(raw_file_path, total_perm_p_below$Group.1, total_perm_p_below$Group.2, total_perm_p_below$empirical_p, NA, NA)
    master_output <- as.data.frame(master_output)

    colnames(master_output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled", "SKAT_O_p-val", "SKAT_O_p-val_resampled")

  }

  master_output
}
