loopSKAT <- function(pre_allocated_SNP_windows, grab_window, window_size = NA, window_shift = NA, window_list = NA,
                     raw_file_path = NA, resampling = FALSE, null_model, n_permutations, parallelize_over = "windows", n_small_null_models = NA,
                     this_scaff_subset, ncore=24, chunk_size=500, backend="foreach", RAM_GB = 4, chunk = TRUE){
  #print("Will use foreach now")




  #iter_object <- icount(length(window_list))

  if(parallelize_over=="windows" & backend=="furrr"){
    #chunk_size <- 1000

    #plan(multiprocess)
    #browser()
    if ( chunk == TRUE ){
      print(paste0("Chunking list of length ",
                   length(pre_allocated_SNP_windows),
                   " into blocks each with no more than",
                   chunk_size,
                   " SNP windows each"))
      pre_allocated_SNP_windows <- split(pre_allocated_SNP_windows, (seq_along(pre_allocated_SNP_windows) - 1) %/% chunk_size)
      print(paste0("...resulting in ", length(pre_allocated_SNP_windows), " blocks"))

    }
    #browser()
    #names(as.list(.GlobalEnv))
    #
    time_to_run_mapping <- proc.time()
    # master_output <- future_map_dfr(.x = pre_allocated_SNP_windows,
    #                                 .f = mappable_SKAT,
    master_output <- future.apply::future_lapply(X = pre_allocated_SNP_windows,
                                                 FUN = mappable_SKAT,
                                                 #.options = future_options(packages = "SKAT"),
                                                 #.progress = TRUE,
                                                 #.cleanup = TRUE,
                                                 #this_position = .x,
                                                 #window_size = window_size,
                                                 #this_scaff_subset = this_scaff_subset,
                                                 raw_file_path = raw_file_path,
                                                 null_model = null_model,
                                                 resampling = resampling,
                                                 n_permutations = n_permutations,
                                                 chunk = chunk)
    message(paste0("Finished parallel run in ",
                 (proc.time() - time_to_run_mapping)[3],
                 "s"))

    master_output_df <- dplyr::bind_rows(master_output)

    #stop()
    #browser()
  }

  if(parallelize_over=="windows" & backend=="foreach"){
    #print("Will parallelize over windows")

    # starting_position <- max((window_size/2),
    #                          round_any(min(this_scaff_subset$POS),
    #                                    window_size/2))
    #
    # window_list <- seq(starting_position,
    #                    max_window,
    #                    by=window_shift)

    master_output <- foreach( i = iter_object,
                              #pos_and_SNPs=ichunk(grab_window, chunkSize=chunkSize),
                              .combine="rbind",
                              .export=c("window_size",
                                        "n_permutations",
                                        "null_model",
                                        "raw_file_path",
                                        "resampling"),
                              .errorhandling="pass",
                              .verbose=TRUE) %do% {

                                pos_and_SNPs = extract_window(this_position = window_list[i],
                                                              window_size = window_size,
                                                              this_scaff_subset = this_scaff_subset)

                                ##print("Inside foreach loop")
                                ##print("Now inside foreach loop. Z is: ")
                                ###print(Z)
                                #stop()
                                ###print("Got Z")

                                if(is.list(pos_and_SNPs) == FALSE){
                                  ##print("No SNPs in this window")
                                  return(rep(-9, 6))
                                }

                                # pos_and_SNPs[[i]]$position
                                # pos_and_SNPs[[i]]$Z

                                # The above function returns None if there are no SNPs within the window. Because of this,
                                # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
                                chunked_output <- as.data.frame(matrix(NA, nrow = 1, ncol = 6))
                                # chunked_output <- foreach(i=1:length(pos_and_SNPs),
                                #                           .combine="rbind"){
                                for(i in 1:length(pos_and_SNPs)){
                                  if(is.matrix(pos_and_SNPs[[i]]$Z)==TRUE){
                                    ##print("About to go into SKAT_one_window")
                                    this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNPs[[i]]$position,
                                                                     window_size,
                                                                     Z = pos_and_SNPs[[i]]$Z,
                                                                     raw_file_path = raw_file_path,
                                                                     resampling = resampling,
                                                                     null_model = null_model,
                                                                     n_permutations = n_permutations)
                                    ##print("Out of SKAT_one_window")
                                  } else {
                                    print(paste0("This chunk's window ", i, " is NOT matrix!"))
                                    this_SKAT_out <- rep(NA, 6)
                                    #browser()
                                  }
                                  chunked_output <- rbind(chunked_output, this_SKAT_out)

                                }
                                ##print("length of SKAT out:")
                                ##print(length(this_SKAT_out))
                                ##print("This SKAT out about to be returned:")
                                ##print(this_SKAT_out)
                                return(chunked_output)
                              }
  }

  if(parallelize_over=="null_models"){
    #print(paste0("Parallelize over ", n_small_null_models, " null models each with n permutations of ", n_permutations, " for ", length(window_list), " windows"))
    master_output <- foreach(null_model_index = 1:n_small_null_models,
                             .combine="rbind",
                             .export=c("window_size",
                                       "n_permutations"),
                             .errorhandling="pass") %dopar% {

                               # if(null_model_index==100){
                               #   #print("Break when index is 100")
                               #   browser()
                               # }

                               #print("Inside foreach loop")
                               #print(paste0("Now inside foreach loop for parallelizing over null models. This # null model is: ", null_model_index))
                               ##print(Z)
                               #stop()
                               ###print("Got Z")

                               #browser()
                               set.seed(null_model_index)
                               null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                                             n.Resampling=n_permutations, type.Resampling="bootstrap")
                               ##print("Made null model")
                               ##browser()
                               master_output <- data.table(matrix(NA, nrow=1, ncol=(n_permutations+2)))
                               #brow


                               ##print("Defined master_output and grab_window on this (to-be) node")
                               ##print(paste0("Length of window_list is: ", length(window_list)))
                               # if(length(window_list) == 6 ){
                               #   browser()
                               # }
                               times <- 0
                               grab_window <- iter_window(this_scaff_subset = this_scaff_subset,
                                                          times = times,
                                                          window_list = window_list,
                                                          window_size = window_size)
                               for(i in 1:length(window_list)){
                                 #print(paste0("On ", i, "th iteration through window_list"))
                                 #pos_and_SNPs <- extract_window(window_list[i], window_size = window_size, this_scaff_subset = this_scaff_subset)
                                 pos_and_SNPs <- nextElem(grab_window)

                                 if(is.list(pos_and_SNPs) == FALSE){
                                   #print("No SNPs in this window")
                                   return(rep(NA, 4))
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
                                   return(rep(NA, 4))
                                 }
                                 #print("length of SKAT out:")
                                 #print(length(this_SKAT_out))
                                 #print("This SKAT out about to be returned:")
                                 #print(this_SKAT_out)
                                 ##browser()
                                 master_output <- rbind(master_output, t(this_SKAT_out))
                                 #print(paste0("Finishing ", i, "th iteration through window_list"))
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
                                 ##print(paste0("For SNP at position ", master_output[1,k]))
                                 n_perm_p_below <- length(which(master_output[ 3:nrow(master_output) , k] <= master_output[ 2 , k ]))
                                 n_perm_p_above <- length(which(master_output[ 3:nrow(master_output) , k] > master_output[ 2 , k ]))
                                 abbreviated_master_output$n_perm_p_below[k] <- n_perm_p_below
                                 abbreviated_master_output$n_perm_p_above[k] <- n_perm_p_above
                               }

                               # if(null_model_index == n_small_null_models){
                               #   browser()
                               # }
                               #print(abbreviated_master_output)
                               if(ncol(abbreviated_master_output) != 4){
                                 #print(paste0("Error: There should be 4 columns for abbreviated master output but instead there are ", ncol(abbreviated_master_output)))
                                 browser()
                                 stop(paste0("Error: There should be 4 columns for abbreviated master output but instead there are ",
                                             ncol(abbreviated_master_output)))
                               }
                               return(abbreviated_master_output)

                             }
    #browser()
    # At this point, we have a data structure containing a row for results from each X permutations results
    # With multiple rows for each SNP. Want to combine rows by SNP and calculate final empirical p-values.

    #master_output$n_perm_p_above <- NULL
    if(nrow(master_output) == 0){
      #print("?????")
      browser()
    }
    #browser()
    total_perm_p_below <- aggregate(master_output$n_perm_p_below,
                                    by=list(master_output$Position, master_output$p_no_resampling), FUN=sum) # https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
    total_perm_p_below$x <- total_perm_p_below$x + 1 # Add one for our p-value without resampling
    total_perm_p_below$empirical_p <- total_perm_p_below$x / ( n_permutations * n_small_null_models )

    master_output <- cbind(raw_file_path, total_perm_p_below$Group.1, total_perm_p_below$Group.2, total_perm_p_below$empirical_p, NA, NA)
    master_output <- as.data.frame(master_output)

    colnames(master_output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled", "SKAT_O_p-val", "SKAT_O_p-val_resampled")

  }
  #browser()

  #stopCluster(cl)
  master_output
}
