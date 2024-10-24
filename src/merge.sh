#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --account=bgmp
#SBATCH --job-name=merge_reads
#SBATCH --output=LOG/merge_reads_%j.out
#SBATCH --error=LOG/merge_reads_%j.err
#SBATCH --time=1-00:00:00            
#SBATCH --nodes=5                  
#SBATCH --cpus-per-task=5
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=ghach@uoregon.edu

conda activate bbmerge
# I really feel like this should be condensed into a script run multiple times with command line arguments
# but I just want to get it working first

set -ue                               # stop on error, if unset variable is accessed
# trailing forward slash on these is not optional, and IDK why
source_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/BLUE/NovaSeq_GC3F_7125/"
source_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/BGMP_2024/RED/NovaSeq_GC3F_7124/"
# This is the lowest directory I can write to in shared/upload
destination_dir1="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/blue/" 
destination_dir2="/projects/bgmp/shared/groups/2024/novel-fluor/shared/upload/NovaSeq_merged/red/" 
# Read in file names in pairs
the_files=("$source_dir1"*.fastq.gz)

mkdir -p "$destination_dir1" 
for ((i=0; i<${#the_files[@]}; i+=2)); do
    # Access the current element and the next one
    if [[ ${the_files[i]} =~ "R1" && ${the_files[i+1]} =~  "R2" ]]; then
        # string slicing to name output files
        sliced_filename=${the_files[i]%%R1*}  # slices string and keeps what occurs before r1
        full_filename="${destination_dir1}${sliced_filename##*NovaSeq_GC3F_7125/}"
        #echo output_filename 
        bbmerge.sh -in1=${the_files[i]} -in2=${the_files[i+1]} -out="${full_filename}MERGED" \
        -outu1="${full_filename}R1_rejected" -outu2="${full_filename}R2_rejected" 
        #echo $full_filename
    else
        echo "Item 1: ${the_files[i]} and Item 2: ${the_files[i+1]} can't be processed."
    fi
done

echo "One down, one to go"

the_other_files=("$source_dir2"*.fastq.gz)

mkdir -p "$destination_dir2" 
for ((i=0; i<${#the_other_files[@]}; i+=2)); do
    # Access the current element and the next one
    if [[ ${the_other_files[i]} =~ "R1" && ${the_other_files[i+1]} =~  "R2" ]]; then
        sliced_filename=${the_other_files[i]%%R1*}  # slices string and keeps what occurs before r1
        full_filename="${destination_dir2}${sliced_filename##*NovaSeq_GC3F_7124/}"
        bbmerge.sh -in1=${the_other_files[i]} -in2=${the_other_files[i+1]} -out="${full_filename}MERGED" \
        -outu1="${full_filename}R1_rejected" -outu2="${full_filename}R2_rejected"
    else
        echo "Item 1: ${the_other_files[i]} and Item 2: ${the_other_files[i+1]} can't be processed."
    fi
done