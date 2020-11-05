calculate_SKAT_empirical_p <- function(Z, n_permutations, null_model,
                                       return_all_p = FALSE, ...){
  this_SKAT_out <- SKAT::SKAT(Z,
                              null_model,
                              # This is the default argument, need not be
                              # explicitly declared
                              #kernel = "linear.weighted",
                              # We wish to allow the user flexibility in
                              #selecting kernel and other options, so pass
                              # arguments with `...`
                              ...)
  if(return_all_p == FALSE){
    p_resampled_SKAT <- ((length(
      subset(
        this_SKAT_out$p.value.resampling,
        this_SKAT_out$p.value.resampling <= this_SKAT_out$p.value)))+1) /
      (n_permutations+1)

    browser()

    return(p_resampled_SKAT)
  }
  if(return_all_p == TRUE){
    #browser()
    return(this_SKAT_out$p.value.resampling)
  }
}
