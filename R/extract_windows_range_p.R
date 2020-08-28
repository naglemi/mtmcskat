extract_windows_range_p <- function(data_in,
                                    upper_bound,
                                    lower_bound){

  # This line should not be needed but is
  # data_in$`SKAT_p-val` <- as.numeric(as.character(data_in$`SKAT_p-val`))

  window_list <- data_in[which(data_in$`SKAT_p-val` < upper_bound & data_in$`SKAT_p-val` >= lower_bound), ]$position
  window_list
}
