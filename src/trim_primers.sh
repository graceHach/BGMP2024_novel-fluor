#!/bin/bash

#SBATCH --partition=<SPECIFY PARTITION>
#SBATCH --account=<SPECIFY PARTITION>
#SBATCH --job-name=trim_primers
#SBATCH --output=LOG/trim_primers_%j.out
#SBATCH --error=LOG/trim_primers_%j.err
#SBATCH --time=10-00:00:00           
#SBATCH --mem=64G

# set up this environment according to setup directions
conda activate ../countbins/

source_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/blue/" 
source_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/red/" 
destination_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/" 
destination_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/red_illum/" 
forward_primers="../primers/CVR205stub_FWD.fasta"
reverse_primers="../primers/CVR205stub_REV.fasta"

# Read in names of merged files
the_files=("${source_dir1}"*_MERGED.fastq)

for file in "${the_files[@]}"; do
    sliced_filename=${file##${source_dir1}}
    echo "processing ${file}"
    hts_Primers -U ${file} -f "${destination_dir1}${sliced_filename%%_MERGED*}" \
   -P ${forward_primers} -Q ${reverse_primers} -l 5 -x -e 6 -d 6 -L "../reports/hts_primers_output/${sliced_filename}" -F
done

echo "Blue complete."

the_files=("${source_dir2}"*_MERGED.fastq)

for file in "${the_files[@]}"; do
    sliced_filename=${file##${source_dir2}}
    echo "processing ${file}"
    hts_Primers -U ${file} -f "${destination_dir2}${sliced_filename%%_MERGED*}" \
   -P ${forward_primers} -Q ${reverse_primers} -l 5 -x -e 6 -d 6 -L "../reports/hts_primers_output/${sliced_filename}" -F
done