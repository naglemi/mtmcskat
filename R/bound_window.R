bound_window <- function(this_position, window_size){
  window_start <- as.numeric(as.character(this_position)) - (window_size/2)

  if (window_start < 0){
    window_start <- 0
  }
  window_end <- as.numeric(as.character(this_position)) + (window_size/2)
  return(list(window_start = window_start, window_end = window_end))
}
