#!/bin/bash


#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=generate_protein_counts
#SBATCH --output=LOG/generate_protein_counts_%j.out
#SBATCH --error=LOG/generate_protein_counts_%j.err     
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

set -eu

in_tsv="../../../shared/dat/counts_per_bin/blue_counts_FINAL.tsv"

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
 -o "../../../shared/dat/counts_per_bin/blue_prot_counts.tsv" -p "../../../shared/dat/NF_pacbio_output/blu/09_final_output/"



echo "starting red"
in_tsv="../../../shared/dat/counts_per_bin/red_counts_FINAL.tsv"

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
 -o "../../../shared/dat/counts_per_bin/red_prot_counts.tsv" -p "../../../shared/dat/NF_pacbio_output/red/09_final_output/"