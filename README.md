
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MTMC-SKAT

<!-- badges: start -->

<!-- badges: end -->

MTMC-SKAT provides high-level R and command line interfaces for SKAT,
most notably with support for multi-threaded adaptive resampling to
rapidly calculate empirical p-values. Resampling is performed using an
algorithm that allows calculation of empirical p-values to a
user-specified number of significant figures. The efficiency of
parallelization over SNP windows or permutations depends on a given
dataset and resources of a given computer. In the automated workflow,
permutations are continually generated and batches of SNP windows are
tested to calculate empirical p-values as sufficient permutations become
available for each batch. The appropriate mode (multithreading over
windows or permutations) is automatically detected for each batch of SNP
windows. See the MTMC-SKAT Workflows vignette for details on this
automatic workflow and information for designing custom workflows.

Workflows are applied over scaffolds independently, with parallelization
over scaffolds or sub-scaffolds on high-performance clusters as
described later. Each scaffold is submitted as a .traw file, which can
be prepared from more common SNP data formats using PLINK. The standard
workflow can be run over a given scaffold using the following command.

``` r
#summary(cars)
```

This function is also accesible from the command line, thus providing a
one-liner that can be easily submitted as a job to any batch query
system.

Finally, the following function can be used to generate a bash script to
submit an array of jobs (each for a given scaffold or sub-scaffold) to
be submitted to to the SGE or SLURMS batch query system.

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
