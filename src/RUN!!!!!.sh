#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=bbmerge_only_parallel_nf
#SBATCH --output=LOG/bbmerge_only_parallel_nf_%j.out
#SBATCH --error=LOG/bbmerge_only_parallel_nf_%j.err        
#SBATCH --nodes=8             
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

# THIS IS THE WHILE SHEBANG. THIS RUNS THE NEXTFLOW PIPELINE THAT PROCESSES ALL ILLUMINA DATA
set -eu
conda activate nextflow-env

/usr/bin/time nextflow run parallel.nf --max_cpus 8

