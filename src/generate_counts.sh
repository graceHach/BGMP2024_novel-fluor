#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=generate_counts
#SBATCH --output=LOG/generate_counts_%j.out
#SBATCH --error=LOG/generate_counts_%j.err     
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

# location of merged files with wildcard expansion 
files="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/*fastq.gz"
# both temp files will be deleted 
intermediate_file="temp.txt"
intermediate_file_2="temp2.txt"

for file in ${files}; do
    echo $file
    # Does this sed expression work on different OSs?
    zcat $file | sed -n '2~4p' >> $intermediate_file
done

# sort and filter for uniqueness
cat $intermediate_file | sort | uniq > $intermediate_file_2
cat $intermediate_file_2 > $intermediate_file
rm $intermediate_file_2

# python script reads in intermediate file to initialize dict
./generate_counts.py --in_dir $files --temp_file $intermediate_file --out_file "../reports/all_blue_counts.tsv"
# delete temp file
rm $intermediate_file