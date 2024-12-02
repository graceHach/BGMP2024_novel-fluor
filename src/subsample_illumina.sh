#!/bin/bash

set -eu
# Take input dir
in_dir_blue="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/BLUE/NovaSeq_GC3F_7125/"
in_dir_red="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/RED/NovaSeq_GC3F_7124/"
out_dir_blue="../mini_illum/blue/"
out_dir_red="../mini_illum/red/"


mkdir -p $out_dir_blue
#BLUE
for file in ${in_dir_blue}*.fastq.gz;
do
    trimmed_filename=(${file##${in_dir_blue}})
    #echo "$trimmed_filename"
    #echo $trimmed_filename
    zcat $file | head -n 10000 > ${out_dir_blue}$trimmed_filename
done

mkdir -p $out_dir_red
#BLUE
for file in ${in_dir_red}*.fastq.gz;
do
    trimmed_filename=(${file##${in_dir_red}})
    #echo "$trimmed_filename"
    #echo $trimmed_filename
    zcat $file | head -n 10000 > ${out_dir_red}$trimmed_filename
done