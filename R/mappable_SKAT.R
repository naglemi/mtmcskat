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

    result_df <- as.data.frame(matrix(NA, nrow=1, ncol=6))

    for ( i in 1:length(pos_and_SNPs)){

      to_append <- SKAT_one_window(this_position = pos_and_SNPs[[i]][[1]],
                                   #window_size,
                                   Z = pos_and_SNPs[[i]][[2]],
                                   raw_file_path = scaffold_ID,
                                   resampling = resampling,
                                   null_model = null_model,
                                   n_permutations = n_permutations)

      result_df <- rbind(result_df, to_append)

    }
    return(result_df[-1,])
  }

  if ( chunk == FALSE ) {

    to_append <- SKAT_one_window(this_position = pos_and_SNPs[[i]][[1]],
                                 #window_size,
                                 Z = pos_and_SNPs[[i]][[2]],
                                 raw_file_path = scaffold_ID,
                                 resampling = resampling,
                                 null_model = null_model,
                                 n_permutations = n_permutations)

    pos_and_SNPs <-NULL
    gc()
    as.data.frame(t(to_append))
  }
}
