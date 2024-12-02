#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=parallel_nf_FINAL
#SBATCH --output=LOG/parallel_nf_FINAL_%j.out
#SBATCH --error=LOG/parallel_nf_FINAL_%j.err        
#SBATCH --nodes=8  
#SBATCH --cpus-per-task=8           
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu
#SBATCH --mem=100G

# THIS IS THE WHILE SHEBANG. THIS RUNS THE NEXTFLOW PIPELINE THAT PROCESSES ALL ILLUMINA DATA
set -eu
conda activate nextflow-env

# COrrect so cores aren't 
nextflow run parallel.nf --max_cpus 8 -resume
#nextflow run parallel.nf --out_dir_trim_red ../test_fastq/similar_filenames
