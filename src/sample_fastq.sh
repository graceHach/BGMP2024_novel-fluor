#!/bin/bash

# This script is intended to be used with time_hts_primers.sh
# Makes a bunch of small fastq files from one large one
destination_dir="../test_fastq/linearly_spaced/"
long_file="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/blue/blue-sort_bin1_S1_L002_MERGED.fastq"

sizes=$(seq 10000000 20000000 160000000)

for size in ${sizes[@]}; do
    # leading zero for lexicographical ordering
    if [ ${size} -lt 100000000 ]
    then
        outfile="${destination_dir}0${size}_lines.fastq"
    else
        outfile="${destination_dir}${size}_lines.fastq"
    fi
    echo ${outfile}
    cat ${long_file} | head -n ${size} > ${outfile}
done