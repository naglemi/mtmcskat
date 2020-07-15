
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Multi-Threaded Monte Carlo Sequence Kernel Association Test (MTMC-SKAT)

<!-- badges: start -->

<!-- badges: end -->

MTMC-SKAT provides high-level R and command line interfaces for
[SKAT](https://github.com/leeshawn/SKAT) (Wu et al., 2011; Lee et al.,
2016), most notably with support for multi-threaded adaptive resampling
to rapidly calculate empirical p-values, and with capabilities for
efficient parallelization both over cores and CPUs. These features can
facilitate high-powered and accurate GWAS for traits that do not satisfy
linear model assumptions.

Resampling is performed using an algorithm that allows calculation of
empirical p-values to a user-specified number of significant figures.
The efficiency of parallelization over SNP windows or permutations
depends on a given dataset and resources of a given computer. In the
automated workflow, permutations are continually generated and batches
of SNP windows are tested to calculate empirical p-values as sufficient
permutations become available for each batch. The appropriate mode
(multithreading over windows or permutations) is automatically detected
for each batch of SNP windows. See the MTMC-SKAT Workflows vignette for
details on this automatic workflow and information for designing custom
workflows.

Workflows are applied over scaffolds independently, with parallelization
over scaffolds or sub-scaffolds on high-performance clusters as
described later. Each scaffold is submitted as a .traw file, which can
be prepared from more common SNP data formats using PLINK. The standard
workflow can be run over a given scaffold using the following command.

``` r
runSKATtraw(phenodata = "poplar_shoot_sample.csv",
            covariates = "poplar_PCs_covariates.csv",
            raw_file_path = "poplar_Chr10_portion.traw",
            window_size = 3000,
            window_shift = 1000,
            output_dir = "Results/",
            job_id = "my_sample_analysis",
            ncore = "AllCores",
            max_accuracy = 5,
            sig_figs = 2)
```

Inputs include file paths to phenotype, covariate and genotype file
paths (in standard PLINK formats for phenotypes and covariates and
`.traw` for genotypes data), the sizes of overlapping SNP windows to be
tested and the shift in position between windows, and the limit for
accuracy (recommended to be set beyond FDR or Bonferroni correction
threshold) given as 10^-x.

This function is also accesible from the command line, thus providing a
one-liner that can be easily submitted as a job to any batch query
system.

``` r
Rscript runSKATtraw.R --phenodata="poplar_shoot_sample.csv" \
covariates="poplar_PCs_covariates.csv" \
raw_file_path="poplar_Chr10_portion.traw" \
window_size=3000 \
window_shift=1000 \
output_dir="Results/" \
job_id="my_sample_analysis" \
ncore="AllCores" \
max_accuracy=5 \
sig_figs=2
```

Finally, the following function can be used to generate a bash script to
submit an array of jobs (each for a given scaffold or sub-scaffold) to
be submitted to to the SGE or SLURMS batch query system…

…

Wu, M.C., Lee, S., Cai, T., Li, Y., Boehnke, M. and Lin, X., 2011.
Rare-variant association testing for sequencing data with the sequence
kernel association test. The American Journal of Human Genetics, 89(1),
pp.82-93.

Lee, S., Fuchsberger, C., Kim, S. and Scott, L., 2016. An efficient
resampling method for calibrating single and gene-based rare variant
association analysis in case–control studies. Biostatistics, 17(1),
pp.1-15.