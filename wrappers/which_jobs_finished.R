library(stringr)
library(tidyr)
library(data.table)

# Which MTMCSKAT jobs finished? Which have csv results; which only have log?
## Want to know which scaffs are too big for memory

results <- list.files("~/Downloads/Results/", recursive = TRUE, pattern = "Pt")
results <- results[!grepl("organized", results)]
# At first, we tried running Chr 5 in 2 parts, before switching to 2
results <- results[!grepl("r5Pt1of2", results)]
results <- results[!grepl("r5Pt2of2", results)]

# Remove some other extra scaffolds
results <- results[!grepl("r2Pt1of2", results)]

results <- results[!grepl("r10Pt1of4", results)]
results <- results[!grepl("r10Pt2of4", results)]
results <- results[!grepl("r10Pt3of4", results)]
results <- results[!grepl("r10Pt4of4", results)]

results <- results[!grepl("r1Pt1of4", results)]
results <- results[!grepl("r1Pt2of4", results)]
results <- results[!grepl("r1Pt3of4", results)]
results <- results[!grepl("r1Pt4of4", results)]

results <- results[!grepl("r4Pt1of2", results)]
results <- results[!grepl("r4Pt2of2", results)]

results <- results[!grepl("r7Pt1of4", results)]
results <- results[!grepl("r7Pt2of4", results)]
results <- results[!grepl("r7Pt3of4", results)]
results <- results[!grepl("r7Pt4of4", results)]

results_log <- results[grep("log", results)]
results_csv <- results[grep("csv", results)]

results_log <- as.data.frame(results_log)
results_csv <- as.data.frame(results_csv)

colnames(results_log)[1] <- "log_path"
colnames(results_csv)[1] <- "csv_path"

results_log$prefix <- tools::file_path_sans_ext(results_log$log_path)
results_csv$prefix <- tools::file_path_sans_ext(results_csv$csv_path)

results <- merge(results_log, results_csv, by = "prefix", all = TRUE)

results$scaff <- stringr::str_split_fixed(results$prefix, "Chr", 2)[, 2]
results$scaff <- stringr::str_split_fixed(results$scaff, "-", 2)[, 1]

data.table::fwrite(results, "~/Downloads/Results/summary_tablenew4.csv")

# Which traits did we run with both P and Q covariates over all chromosomes?

traits <- str_split_fixed(results[, 3], "-", 3)[, 2]

results$traits <- traits

traits <- unique(traits)
traits <- traits[traits != ""] # Get rid of blank that wound up in there

results <- results[results$traits != "", ]

results_table <- as.data.frame(table(results$traits, results$scaff))

colnames(results_table)[1:2] <- c("trait", "scaffold")

table(results$traits)

results_table_wide <- spread(results_table,
                             scaffold,
                             Freq)

fwrite(results_table_wide, "~/Downloads/Results/results_table_wide.csv")

results_table_all_scaff <- as.data.frame(table(results$traits))

traits_studied_both_P_and_Q <-
  results_table_all_scaff$Var1[which(results_table_all_scaff$Freq == 104)]


# Reorganize files --------------------------------------------------------

getwd()

if (!dir.exists("~/Downloads/Results/organized/")){
  dir.create("~/Downloads/Results/organized/")
}

if (!dir.exists("~/Downloads/Results/organized/SLURMS_7K")){
  dir.create("~/Downloads/Results/organized/SLURMS_7K")
}

if (!dir.exists("~/Downloads/Results/organized/SLURMS_6PC")){
  dir.create("~/Downloads/Results/organized/SLURMS_6PC")
}

setwd("~/Downloads/Results/")

for(trait in traits_studied_both_P_and_Q){
  print(trait)
  csvs_this_trait <- results$csv_path[grepl(trait, results$csv_path)]
  print(length(csvs_this_trait))
  for(csv in csvs_this_trait){
    if(grepl("7K", csv)){
      file.copy(csv, paste0("~/Downloads/Results/organized/SLURMS_7K/",
                            basename(csv)))
    }
    if(grepl("6PC", csv)){
      file.copy(csv, paste0("~/Downloads/Results/organized/SLURMS_6PC/",
                            basename(csv)))
    }
  }
}

# Which traits did we run at least P covariates over all chromosom --------

results_p_only <- results[grepl("6PC", results$prefix), ]

results_table <- as.data.frame(table(results_p_only$traits,
                                     results_p_only$scaff))

colnames(results_table)[1:2] <- c("trait", "scaffold")

table(results_p_only$traits)

results_table_wide <- spread(results_table,
                             scaffold,
                             Freq)

fwrite(results_table_wide, "~/Downloads/Results/results_p_only_table_widea2.csv")

results_table_all_scaff <- as.data.frame(table(results_p_only$traits))

traits_studied_P <-
  results_table_all_scaff$Var1[which(results_table_all_scaff$Freq == 52)]


# More file reorganization, for those studied with P only -----------------

for(trait in traits_studied_P){
  print(trait)
  csvs_this_trait <- results_p_only$csv_path[grepl(trait,
                                                   results_p_only$csv_path)]
  print(length(csvs_this_trait))
  for(csv in csvs_this_trait){
    # if(grepl("7K", csv)){
    #   file.copy(csv, paste0("~/Downloads/Results/organized/SLURMS_7K/",
    #                         basename(csv)))
    # }
    if(grepl("6PC", csv)){
      file.copy(csv, paste0("~/Downloads/Results/organized/SLURMS_6PC/",
                            basename(csv)))
    }
  }
}

