mappable_SKAT <- function(#this_position,
                          #these_positions,
                          pos_and_SNPs,
                          window_size,
                          this_scaff_subset,
                          raw_file_path,
                          null_model,
                          n_permutations,
                          resampling=FALSE,
                          return_all_p_vals=FALSE,
                          chunk = TRUE){

  if ( chunk == TRUE ) {

    #print(paste0("from inside mappable_SKAT, ", length(pos_and_SNPs), "positions from ", pos_and_SNPs[[1]][1], " to ", pos_and_SNPs[[length(pos_and_SNPs)]][1]))

    result_df <- as.data.frame(matrix(NA, nrow=1, ncol=6))
    #browser()
    for ( i in 1:length(pos_and_SNPs)){
      #print(paste0("In this chunk (list of positions and Z matrices) we are at index ", i, " of ", length(pos_and_SNPs)))
      #browser()      print(paste0("Dim of "))
      # print(paste0("SNP sublist has length of ", length(pos_and_SNPs)))
      # print(dim(pos_and_SNPs[[i]][[2]]))
      # print("Is it a matrix?")
      # print(is.matrix(pos_and_SNPs[[i]][[2]]))
      # print(paste0(" at ", pos_and_SNPs[[i]][[1]]))
      to_append <- callSKAT(pos_and_SNPs = pos_and_SNPs[[i]], scaffold_ID = raw_file_path, resampling = TRUE, null_model = null_model, n_permutations = n_permutations)
      # print(paste0("Raw file path is: ", raw_file_path))
      # print("Will append....")
      # #browser()
      # print(to_append)
      result_df <- rbind(result_df, to_append)

    }
    return(result_df[-1,])
  }

  if ( chunk == FALSE ) {
    #print("Chunking is off. This is deprecated and might not work any more. Functionality to be removed.")

    to_append <- callSKAT(pos_and_SNPs = pos_and_SNPs, scaffold_ID = raw_file_path, resampling = TRUE, null_model = null_model, n_permutations = n_permutations)

    #print(length(to_append))
    ##print(to_append)
    pos_and_SNPs <-NULL
    gc()
    as.data.frame(t(to_append))
  }
}
