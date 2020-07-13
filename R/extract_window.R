extract_window <- function(this_position, window_size, this_scaff_subset){
  #print("Just entered extract_window")
  #print(paste0("This position (within extract_window) is: ", this_position, " with SNPs: "))

  if(is.na(this_position)){
    #print("NA position. Reached end of scaffold?")
    #print("Passing StopIteration")
    #stop("StopIteration")
    return(this_position, NA)
  }

  locus_of_interest <- this_position
  window_start <- as.numeric(as.character(locus_of_interest)) - (window_size/2)

  if (window_start < 0){
    window_start <- 0
  }
  window_end <- as.numeric(as.character(locus_of_interest)) + (window_size/2)

  #print(paste0("SNP window: ", window_start, "-", window_end))

  indices_to_pull <- which(with(this_scaff_subset, this_scaff_subset$POS >= window_start & this_scaff_subset$POS <= window_end), arr.ind = TRUE)
  #print("Got indices if they exist")
  if (length(indices_to_pull) == 0){
    #this_position <- this_position + window_shift
    #next
    #return()
    # I think we need return(NA) instead of return() to stop is.matrix(Z) from giving a "missing value where TRUE/FALSE needed" error
    return(list(this_position, NA))
  }
  #print("Determined whether there are any indices")
  genodata_thiswindow <- this_scaff_subset[indices_to_pull[1]:indices_to_pull[length(indices_to_pull)],]
  #return(genodata_thiswindow)
  #print("Obtained subseted genodata")
  #browser()

  Chr <- unique(genodata_thiswindow[,1])
  if(length(Chr) < 1) {
    stop("Where is chromosome data?")
  }
  if(length(Chr) > 1) {
    stop(paste0(Sys.time(),
                " - extract_window should be provided with no more than a single scaffold, but appears to have multiple."))
  }

  # genodata_thiswindow[,1:6] <- NULL
  # genodata_thiswindow <- data.frame(genodata_thiswindow)
  # Z <- t(as.matrix(genodata_thiswindow))
  # colnames(Z) <- NULL
  Z <- convert_to_Z(genodata_thiswindow = genodata_thiswindow)

  #print("Z extracted dimensions and first 10 col, row:")
  #print(dim(Z))
  #print(head(Z)[,1:min(10, ncol(Z))])
  #print("Now about to return Z")
  #browser()
  #browser()


  return(list(this_position, Z, Chr))
}
