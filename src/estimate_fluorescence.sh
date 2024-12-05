#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=estimate_fluorescence
#SBATCH --output=LOG/estimate_fluorescence_%j.out
#SBATCH --error=LOG/estimate_fluorescence_%j.err     
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

# This script replaces generate_counts_better.sh, generate_protein_counts.sh, 
# counts_matrix_manipulation.sh, and bins_to_MEF.sh

set -eu

#######       GENERATE ILLUMINA COUNTS       #######   
# IS THE WAY TO DO THIS IS NEXTFLOW TO HAVE THE --IN_DIR BELOW BE ITSELF A COMMAND-LINE PARAMETER?  
#./generate_counts_better.py --in_dir "/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/" --out_file "../reports/blue_counts.tsv"
#./generate_counts_better.py --in_dir "/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/red_illum/" --out_file "../reports/red_counts.tsv"

#######       INCORPORATE PROTEIN SEQUENCES       ####### 
in_tsv="../reports/blue_counts.tsv"

column_sums=$(awk '{
    for (i=2; i<=NF; i++) {  
        sum[i] += $i
    }
} END {
    for (i=1; i<=NF; i++) {
        print sum[i]
    }
}' ${in_tsv})

./generate_protein_counts.py -cs $column_sums -f ${in_tsv} \
 -o "../reports/blue_prot_counts.tsv" -p "../../../shared/dat/NF_results/BLUE/09_final_output/"


in_tsv="../reports/red_counts.tsv"

column_sums=$(awk '{
    for (i=2; i<=NF; i++) {  
        sum[i] += $i
    }
} END {
    for (i=1; i<=NF; i++) {
        print sum[i]
    }
}' ${in_tsv})

./generate_protein_counts.py -cs $column_sums -f ${in_tsv} \
 -o "../reports/red_prot_counts.tsv" -p "../../../shared/dat/NF_results/RED/09_final_output/"

#######       DO THE MATH TO GET MEDIAN BINS
in_file="../reports/blue_prot_counts.tsv"
out_file="../reports/blue_median_bin.tsv"

./counts_matrix_manipulation.py --in_tsv $in_file --out_file $out_file

in_file="../reports/red_prot_counts.tsv"
out_file="../reports/red_median_bin.tsv"

./counts_matrix_manipulation.py --in_tsv $in_file --out_file $out_file

#######       BINS TO MEF
./bins_to_MEF.py --in_bin_edges "../bins/blue_bins.tsv"  --in_median_bins "../reports/blue_median_bin.tsv" --out_file "../reports/blue_MEF_FINAL.tsv"

./bins_to_MEF.py --in_bin_edges  "../bins/red_bins.tsv" --in_median_bins "../reports/red_median_bin.tsv" --out_file "../reports/red_MEF_FINAL.tsv"
