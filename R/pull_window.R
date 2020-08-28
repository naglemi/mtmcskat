pull_window <- function(this_scaff_subset, window_start, window_end) {
  indices_to_pull <- which(with(this_scaff_subset, this_scaff_subset$POS >= window_start & this_scaff_subset$POS <= window_end), arr.ind = TRUE)
  #print("Got indices if they exist")
  if (length(indices_to_pull) == 0){
    #this_position <- this_position + window_shift
    #next
    #return()
    # I think we need return(NA) instead of return() to stop is.matrix(Z) from giving a "missing value where TRUE/FALSE needed" error
    return(NULL)
  }
  #print("Determined whether there are any indices")
  genodata_thiswindow <- this_scaff_subset[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)],]
  return(genodata_thiswindow)
}
