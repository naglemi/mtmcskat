map_SKAT_nm <- function(null_model_index,
                        pos_and_SNP_list,
                        window_size,
                        raw_file_path,
                        n_permutations,
                        resampling=TRUE){

  set.seed(null_model_index)
  null_model <- SKAT_Null_Model(this_phenotype ~ 1 + covariates$V1 + covariates$V2 + covariates$V3 + covariates$V4 + covariates$V5 + covariates$V6 + covariates$V7 + covariates$V8 + covariates$V9 + covariates$V10 + covariates$V11, out_type="C",
                                n.Resampling=n_permutations, type.Resampling="bootstrap")
  ##print("Made null model")
  ##browser()
  chunked_output <- data.table(matrix(NA, nrow=1, ncol=(n_permutations+2)))

  #print(pos_and_SNPs[[1]])

  # The above function returns None if there are no SNPs within the window. Because of this,
  # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
  for(i in 1:length(pos_and_SNP_list)){

    if(is.matrix(pos_and_SNP_list[[i]][[2]])==TRUE){
      #print("About to go into SKAT_one_window")
      this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNP_list[[i]][[1]],
                                       window_size,
                                       Z = pos_and_SNP_list[[i]][[2]],
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
    chunked_output <- rbind(chunked_output, t(this_SKAT_out))
    #print(paste0("Finishing ", i, "th iteration through window_list"))
  }

  p_null_tallies <- tally_p_null(chunked_output[-1,])

  return(p_null_tallies)
}
