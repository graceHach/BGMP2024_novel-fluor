#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=generate_counts_better
#SBATCH --output=LOG/generate_counts_better_%j.out
#SBATCH --error=LOG/generate_counts_better_%j.err     
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

set -eu
# python script reads in intermediate file to initialize dict
./generate_counts_better.py --in_dir "/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/" --out_file "../reports/blue_counts_FINAL.tsv"
./generate_counts_better.py --in_dir "/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/red_illum/" --out_file "../reports/red_counts_FINAL.tsv"