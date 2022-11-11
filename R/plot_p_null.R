plot_p_null <- function(Z, null_model, missing_cutoff,
                        output_basename, output_dir,
                        scaffold_ID, this_position){

  #browser()

  this_SKAT_out <- SKAT::SKAT(Z,
                              null_model,
                              missing_cutoff = 0.15
                              # This is the default argument, need not be
                              # explicitly declared
                              #kernel = "linear.weighted",
                              # We wish to allow the user flexibility in
                              #selecting kernel and other options, so pass
                              # arguments with `...`
  )

  if(is.na(output_basename)){
    output_basename <- "Untitled_
        job"
  }
  if(is.na(output_dir)){
    output_dir <- "plot_p_null"
  }
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  file_dest <- paste0(output_dir, "/", # "/", output_basename, "_pnulltoplot_",
                      scaffold_ID, "_", this_position, ".csv")
  #browser()
  print(paste0("Writing p-vals to: ", file_dest, " from ", getwd()))
  if(!file.exists(file_dest)){
    p_to_write <- sample(this_SKAT_out$p.value.resampling,
                         min(length(this_SKAT_out$p.value.resampling), 10000))
    data.table::fwrite(as.data.frame(p_to_write),
                       file_dest, col.names = FALSE)
  }
}
