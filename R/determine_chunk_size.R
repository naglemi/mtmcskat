determine_chunk_size <- function(cores = detectCores(),
                                 n_jobs){
  # We want to each core to get 3-10 jobs, a good range that helps avoid too much master-slave communication
  # and balances the load well

}

determine_switch_point <- function(Null_Model, available_RAM){
  # We want to switch modes once there is no longer enough RAM for a null model that is large enough.
  # This sholud be integrated with a high-level "RAM usage" option that also integrates with option() for setting future global memory allocation... will involve a calculation with the number of available cores, memory per core...
}
