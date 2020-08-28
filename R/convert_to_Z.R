convert_to_Z <- function(genodata_thiswindow){
  if(ncol(genodata_thiswindow) < 6){
    print("Should be at least 6 columns here")
    browser()
  }
  genodata_thiswindow[,1:6] <- NULL
  genodata_thiswindow <- data.frame(genodata_thiswindow)

  Z <- t(as.matrix(genodata_thiswindow))
  colnames(Z) <- NULL
  return(Z)
}
