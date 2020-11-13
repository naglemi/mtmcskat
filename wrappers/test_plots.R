results_path <- "/scratch2/NSF_GWAS/Results/SKAT/"


# Shoot phenotype at week 5 -----------------------------------------------

results_pattern <- "pheno-shoot_5w"

results_files <- list.files(results_path,
                            pattern = results_pattern,
                            full.names = TRUE)

results_files[c(5,7)]

library(data.table)
file1 <- fread(results_files[5])
file2 <- fread(results_files[7])

file_combined <- rbind(file1, file2)

master_output <- file_combined

par(mfrow=c(2,1))

qqman::manhattan(
  master_output[!is.na(master_output$`SKAT_p-val`), ],
  chr = "Chr",
  bp = "position",
  p = "SKAT_p-val",
  ylim = c(0,
           ceiling(max(
             na.omit((-log(master_output$`SKAT_p-val`,
                           base = 10)))))),
  suggestiveline = FALSE,
  genomewideline = FALSE)

qqman::manhattan(
  master_output[!is.na(master_output$`SKAT_p-val_resampled`), ],
  chr = "Chr",
  bp = "position",
  p = "SKAT_p-val_resampled",
  ylim = c(0,
           ceiling(max(
             na.omit((-log(master_output$`SKAT_p-val`,
                           base = 10)))))),
  suggestiveline = FALSE,
  genomewideline = FALSE)


# Callus phenotype at week 5 ----------------------------------------------

results_pattern <- "pheno-callus_5w"

results_files <- list.files(results_path,
                            pattern = results_pattern,
                            full.names = TRUE)

results_files

library(data.table)
file1 <- fread(results_files[1])
file2 <- fread(results_files[3])

file_combined <- rbind(file1, file2)

master_output <- file_combined

par(mfrow=c(2,1))

qqman::manhattan(
  master_output[!is.na(master_output$`SKAT_p-val`), ],
  chr = "Chr",
  bp = "position",
  p = "SKAT_p-val",
  ylim = c(0,
           ceiling(max(
             na.omit((-log(master_output$`SKAT_p-val`,
                           base = 10)))))),
  suggestiveline = FALSE,
  genomewideline = FALSE)

qqman::manhattan(
  master_output[!is.na(master_output$`SKAT_p-val_resampled`), ],
  chr = "Chr",
  bp = "position",
  p = "SKAT_p-val_resampled",
  ylim = c(0,
           ceiling(max(
             na.omit((-log(master_output$`SKAT_p-val`,
                           base = 10)))))),
  suggestiveline = FALSE,
  genomewideline = FALSE)
