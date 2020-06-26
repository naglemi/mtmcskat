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
                                                                null_model = null_model_noresample,
                                                                n_permutations = n_permutations)
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
