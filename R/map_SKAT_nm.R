map_SKAT_nm <- function(null_model_index,
                        this_phenotype,
                        covariates,
                        pos_and_SNP_list,
                        scaffold_ID,
                        n_permutations,
                        resampling=TRUE){

  set.seed(null_model_index)

  null_model <- SKAT::SKAT_Null_Model(
    this_phenotype ~ 1 + as.matrix(covariates), out_type="C",
    n.Resampling=n_permutations,
    type.Resampling="bootstrap")

  chunked_output <- data.table::data.table(
    matrix(NA, nrow=1, ncol=(n_permutations+2)))

  # The above function returns None if there are no SNPs within the window. Because of this,
  # we only run the below if a matrix is returned to avoid "Z is not a matrix" error.
  for(i in 1:length(pos_and_SNP_list)){

    if(is.matrix(pos_and_SNP_list[[i]][[2]])==TRUE){
      #print("About to go into SKAT_one_window")
      this_SKAT_out <- SKAT_one_window(this_position = pos_and_SNP_list[[i]][[1]],
                                       Z = pos_and_SNP_list[[i]][[2]],
                                       scaffold_ID = scaffold_ID,
                                       resampling = resampling,
                                       null_model = null_model,
                                       n_permutations = n_permutations,
                                       return_all_p_vals = TRUE)

    } else { # This condition should never happen since we next if no SNPs
      return(rep(NA, 4))
    }

    chunked_output <- rbind(chunked_output, t(this_SKAT_out))
  }

  p_null_tallies <- tally_p_null(chunked_output[-1,])

  null_model <- NULL
  gc()

  return(p_null_tallies)
}
