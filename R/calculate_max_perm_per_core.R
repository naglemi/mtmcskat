calculate_max_perm_per_core <- function(nm_RAM_per_perm,
                                        RAM,
                                        ncore){
  # What is the maximum number of permutations that can be performed without a round of boss-worker communication?
  max_simultaneous_perm <- as.numeric(RAM / nm_RAM_per_perm)
  #browser()

  max_simultaneous_perm_per_core <- ceiling(max_simultaneous_perm / ncore)

  max_simultaneous_perm_per_core
}



benchmark_nm <- function(phenotype, covariates, benchmark_size = 1000){
  # benchmark_nm(phenotype = this_phenotype, covariates = covariates)
  null_model <- SKAT_Null_Model(phenotype ~ 1 + as.matrix(covariates), out_type="C",
                                n.Resampling = benchmark_size, type.Resampling="bootstrap")
  nm_RAM_per_perm <- object.size(null_model) / benchmark_size
  nm_RAM_per_perm
}



size_RAM_partition <- function(partition_factor){
  # size_RAM_partition(2)
  RAM <- as.numeric(benchmarkme::get_ram())/partition_factor
  RAM
}

