calculate_max_perm_per_core <- function(nm_RAM_per_perm,
                                        RAM,
                                        n_core){
  # What is the maximum number of permutations that can be performed without a round of boss-worker communication?
  max_simultaneous_perm <- as.numeric(RAM / nm_RAM_per_perm)
  #browser()

  max_simultaneous_perm_per_core <- floor(max_simultaneous_perm / n_core)

  max_simultaneous_perm_per_core
}



benchmark_nm <- function(phenotype, covariates, benchmark_size = 1000){
  # benchmark_nm(phenotype = this_phenotype, covariates = covariates)
  null_model <- SKAT::SKAT_Null_Model(phenotype ~ 1 + as.matrix(covariates),
                                      out_type="C",
                                      n.Resampling = benchmark_size,
                                      type.Resampling="bootstrap")

  nm_RAM_per_perm <- object.size(null_model) / benchmark_size
  nm_RAM_per_perm
}

size_RAM_wiggle <- function(wiggle_factor){
  # size_RAM_wiggle(2)
  RAM <- as.numeric(benchmarkme::get_ram())/wiggle_factor
  RAM
}

