p_empirical_from_tally <- function(p_null_tallies, scaffold_ID){
  total_perm_p_ltoreq <- stats::aggregate(
    p_null_tallies$n_perm_ltoreq,
    by=list(p_null_tallies$position,
            p_null_tallies$`SKAT_p-val`),
    FUN=sum) # https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group

  total_perm_p_ltoreq$x <- total_perm_p_ltoreq$x + 1 # Add one for our p-value without resampling

  total_perm_p_above <- stats::aggregate(
    p_null_tallies$n_perm_above,
    by=list(p_null_tallies$position,
            p_null_tallies$`SKAT_p-val`), FUN=sum)

  total_perm_p_ltoreq$empirical_p <- total_perm_p_ltoreq$x /
    ( total_perm_p_ltoreq$x + total_perm_p_above$x )

  output <- cbind(scaffold_ID,
                  total_perm_p_ltoreq$Group.1,
                  total_perm_p_ltoreq$Group.2,
                  total_perm_p_ltoreq$empirical_p)

  output <- as.data.frame(output)

  colnames(output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled")

  output$`SKAT_p-val` <- as.numeric(as.character(output$`SKAT_p-val`))
  output$`SKAT_p-val_resampled` <- as.numeric(as.character(
    output$`SKAT_p-val_resampled`))

  output
}
