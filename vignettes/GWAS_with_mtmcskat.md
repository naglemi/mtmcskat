In this vignette, we describe key mid-level functions provided by the
`MTMC-SKAT` package, their usage in the standard automated workflow and
how they can be used to produce a custom workflow.

    #library(SKATMCMT)

The function `pre_allocate` is used to load SNP data from the `.traw`
format and convert it to a list of SNP windows. The resulting list is
saved in the directory specified by `pre_allocated_dir`. To avoid
regenerating a list used for multiple analyses, `pre_allocate` will read
in the corresponding list for a genotype file if one exists.

    pre_allocated_SNP_windows <- pre_allocate(raw_file_path = "poplar_Chr10_portion.traw",
                                              window_size = 3000,
                                              window_shift = 1000,
                                              pre_allocated_dir = "pre_allocated/")

Upon pre-allocating SNP window data, we can begin the initial round of
SKAT, without resampling, to produce inital p-values for each SNP window
without any permutations. From these initial results, the user or
preprogrammed filters can select SNP windows deemed interesting and
perform `MTMC-SKAT` to calculate empirical p-values for these SNP
windows. The provided function for multi-threaded SKAT (`mtskat`) runs
SKAT over pre-allocated SNP windows, with chunks of SNP windows
distributed among cores.

    results_sans_resampling <- mtskat(pre_allocated_SNP_windows = pre_allocated_SNP_windows
                                      phenotype = read.csv("poplar_shoot_sample.csv"),
                                      covariates = read.csv("poplar_PC_covariates.csv"),
                                      ncore = "AllCores")

Based on these initial results, batches of SNP windows producing
p-values over a given range of interest can be selected using
`extract_windows_range_p`. As in the standard automated workflow, we
recommend the range of p-values for each batch of SNPs vary by no more
than an order of magnitude, to allow calculation of empirical p-values
to a desired accuracy without an unnecesary number of permutations.

    window_list <- extract_windows_range_p(data = results_sans_resampling,
                                           upper_bound = 10^-2,
                                           lower_bound = 10^-3)
    new_pre_allocated_SNP_windows <- subset_list_of_lists(list_of_lists = pre_allocated_SNP_windows,
                                                          desired_list = window_list)

Once this subset of SNP windows has been pre-allocated into a new list,
empirical p-values can be calculated using `mtmcskat` with the
`multithreading` option set to `snp` (to parallelize over SNP windows)
or `nm` (to parallelize over permuted null models). The former is more
efficient when 1) a sufficient number of permuted null models for the
desired accuracy of empirical p-values can be stored in memory; and 2)
The number of SNP windows is no less than the number of cores. For lower
p-value ranges requiring larger number of permutations and with smaller
batches of SNP windows, the latter mode becomes more efficient. In the
standard automated workflow, the appropriate function is used for each
batch of SNPs based on these criteria.

    batch_results_with_resampling_over_SNPs <- mtmcskat(pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
                                                        phenotype = read.csv("poplar_shoot_sample.csv"),
                                                        covariates = read.csv("poplar_PC_covariates.csv"),
                                                        ncore = "AllCores",
                                                        multithreading = "snp")

    batch_results_with_resampling_over_NMs <- mtmcskat(pre_allocated_SNP_windows = new_pre_allocated_SNP_windows,
                                                        phenotype = read.csv("poplar_shoot_sample.csv"),
                                                        covariates = read.csv("poplar_PC_covariates.csv"),
                                                        ncore = "AllCores",
                                                        multithreading = "nm")

In rare cases, empirical p-values may be significantly lower than
initial p-values, leading to fewer significant figures for an empirical
p-value than expected given the number of permutations run. These SNP
windows can then be extracted for subsequent rounds of analyses with
more permutations, until the desired number of significant figures or
accuracy limit is reached.

    window_list <- extract_windows_sig_fig(data = batch_results_with_resampling_over_NMs,
                                           n_sig_fig = 1) # Suppose we desire two significant figures and thus want all windows for which empirical p-values only have one
    new_pre_allocated_SNP_windows <- subset_list_of_lists(list_of_lists = pre_allocated_SNP_windows,
                                                          desired_list = window_list)

…
