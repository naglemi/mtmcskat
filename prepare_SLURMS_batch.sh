#!/bin/bash

while getopts d:p:c:i:h:o:j:n:m:s:t: flag
do
 case "${flag}" in
  d) data_folder=${OPTARG};;
  p) phenodata=${OPTARG};;
  c) covariates=${OPTARG};;
  i) window_size=${OPTARG};;
  h) window_shift=${OPTARG};;
  o) output_dir=${OPTARG};;
  j) job_id=${OPTARG};;
  n) ncore=${OPTARG};;
  m) max_accuracy=${OPTARG};;
  s) sig_figs=${OPTARG};;
  t) max_time=${OPTARG};;
 esac
done

echo "Preparing script for submitting array of jobs to SLURMS scheduler..."
echo "Folder of scaffolds, each for one job: $data_folder"
echo "Phenotype file: $phenodata"
echo "Covariate file: $covariates"
echo "SNP window size: $window_size"
echo "SNP window shift: $window_shift"
echo "Directory to save results to: $output_dir"
echo "Job identifier: $job_id"
echo "Number of threads (per job): $ncore"
echo "Maximum accuracy for empirical p-values: $max_accuracy"
echo "Desired # significant figures for empirical p-values: $sig_figs"
echo ""

job_list_name="${job_id}.sh"

if test -f "$job_list_name"; then
 echo "WARNING: job file $job_list_name already exists... so this job file will be created with the name ${job_id}_$(date '+%H%M%S%m%d%Y').sh"
 echo ""
fi
 
job_list_name="${job_id}_$(date '+%H%M%S%m%d%Y').sh"

max_chr=$(ls $data_folder | wc -l)

echo "Number of scaffolds (equivalent to # jobs): $max_chr"
echo ""

set +H
echo "#!/bin/bash" >> $job_list_name

echo "#SBATCH --job-name=$phenodata" >> $job_list_name
echo "#SBATCH --output=$job_name_%A_%a.out" >> $job_list_name
echo "#SBATCH --error=$job_name%A_%a.err" >> $job_list_name
echo "#SBATCH --array=1-$max_chr" >> $job_list_name
echo "#SBATCH --time=$max_time" >> $job_list_name
echo "module load R"

for file in $data_folder/*
do
scaffold_name=$(basename $file)
jobname="${job_id}_$(basename ${scaffold_name%.*})"
echo "Adding job for input and output:"
echo $file
echo $jobname
echo "Rscript MTMCSKAT_workflow.R --phenodata=$phenodata \
--covariates=$covariates \
--raw_file_path=$file \
--window_size=$window_size \
--window_shift=$window_shift \
--output_dir=$output_dir \
--job_id=$jobname \
--ncore=$ncore \
--max_accuracy=$max_accuracy \
--sig_figs=$sig_figs" >> $job_list_name
done
