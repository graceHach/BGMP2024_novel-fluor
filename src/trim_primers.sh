#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=trim_primers
#SBATCH --output=LOG/trim_primers_%j.out
#SBATCH --error=LOG/trim_primers_%j.err
#SBATCH --time=10-00:00:00           
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu
#SBATCH --mem=64G

conda activate htstream

source_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/blue/" 
source_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/red/" 
destination_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/" 
destination_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/red_illum/" 
forward_primers="../primers/CVR205stub_FWD.fasta"
reverse_primers="../primers/CVR205stub_REV.fasta"

# Read in names of merged files
the_files=("${source_dir1}"*_MERGED.fastq)

echo "I beseech you, in the bowels of Christ, think it possible you may be mistaken."

for file in "${the_files[@]}"; do
    sliced_filename=${file##${source_dir1}}
    echo "processing ${file}"
    hts_Primers -U ${file} -f "${destination_dir1}${sliced_filename%%_MERGED*}" \
   -P ${forward_primers} -Q ${reverse_primers} -l 5 -x -e 6 -d 6 -L "../reports/hts_primers_output_BLUE" -F
done

echo "Blue complete."

the_files=("${source_dir2}"*_MERGED.fastq)

for file in "${the_files[@]}"; do
    sliced_filename=${file##${source_dir2}}
    echo "processing ${file}"
    hts_Primers -U ${file} -f "${destination_dir2}${sliced_filename%%_MERGED*}" \
   -P ${forward_primers} -Q ${reverse_primers} -l 5 -x -e 6 -d 6 -L "../reports/hts_primers_output_RED" -F
done