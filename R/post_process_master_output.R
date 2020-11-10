post_process_master_output <- function(master_output){
  if(is.na(master_output)){
    print("NA output")
    return(NA)
  }
  master_output <- data.frame(master_output)
  colnames(master_output) <- c("Chr", "position", "SKAT_p-val", "SKAT_p-val_resampled")
  master_output$`SKAT_p-val` <- as.numeric(as.character(master_output$`SKAT_p-val`))
  master_output$`SKAT_p-val_resampled` <- as.numeric(as.character(master_output$`SKAT_p-val_resampled`))
  master_output
}
