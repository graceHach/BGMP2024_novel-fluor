#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=time_hts_primers
#SBATCH --output=timing_LOG/time_hts_primers_%j.out
#SBATCH --error=timing_LOG/time_hts_primers_%j.err     
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

file_dir="../test_fastq/linearly_spaced/"
destination_dir1="../test_fastq/linearly_spaced_out/" 
forward_primers="../primers/CVR205stub_FWD.fasta"
reverse_primers="../primers/CVR205stub_REV.fasta"
output_file="../reports/hts_primers_runtime_10_150_million.tsv"

mkdir -p ${destination_dir1}


set -eu
conda activate htstream

echo -e "lines in the fastq\tReal time elapsed (hh:mm:ss)" > ${output_file}

# no parenthesis, $() makes it interpet the filename as a command
# pure glob
files=${file_dir}*.fastq

for file in ${files[@]}; do
    lines=$(wc ${file} -l)
    sliced_filename=${file##${file_dir}}
    wc ${file} -l | awk '{printf "%s\t", $1}' >> ${output_file}
    # time actually directs to standard error. 
    # output=$( { time COMMAND; } 2>&1 ) redirects cobines stdout and stderr into one stream
    # wasn't working w/ just time -f. Printing elapsed time only
    /usr/bin/time -f "%E" hts_Primers -U ${file} -f "${destination_dir1}${sliced_filename}" \
    -P ${forward_primers} -Q ${reverse_primers} -l 5 -x -e 6 -d 6 \
    -F >> ${output_file} 2>&1
    # above directs standard out to output file, THEN merges standard error and standard out
done