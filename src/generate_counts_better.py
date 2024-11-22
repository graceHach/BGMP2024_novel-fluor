#!/usr/bin/env python

import argparse
import gzip
import glob
import numpy as np

# Run using helper script generate_counts.sh

# input directory, output directory 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--in_dir",type=str,help="directory containing merged, primer trimmed reads", default="/projects/bgmp/shared/groups/2024/novel-fluor/shared/dat/blue_illum/")
    parser.add_argument("-o","--out_file",type=str,help="", default="../reports/blue_counts_better_script.tsv")
    return parser.parse_args()

args = get_args()

not_20_bp = 0
# read in intermmediate file, initialize the dictionary
print("initializing counts dict", flush=True)
counts = {} 
# keys are UMIs, values are numpy arrays with nine entries

files = glob.glob(args.in_dir+"*.fastq.gz")
# glob puts these in a weird order 
files_list = sorted(list(files))
# sort the files???
for file_index, file in enumerate(files_list):
    print("processing: ",file, flush=True)
    with gzip.open(file, 'rt') as fh:
        for index, line in enumerate(fh):
            if index%4==1:
                if len(line)==21:
                    if line[:-1] in counts:
                        counts[line[:-1]][file_index] += 1
                    else:
                        # I USED INT8 INITIALLY AND IT FRICKING UNDERFLOWED AND HALF THE NUMBERS WERE NEGATIVE AAAAAUGH
                        # 4.5 HOUR RUNTIME 
                        counts[line[:-1]] = np.zeros(9, dtype="int64")
                        counts[line[:-1]][file_index] = 1
                else:
                    not_20_bp += 1

# write list of file names as header of file
header_line = "Barcode\t" + '\t'.join(files_list) + "\n"

print("Starting writing")
keys = counts.keys()
with open(args.out_file, 'wt') as fh_out:
    fh_out.write(header_line)
    for key in keys:
        write_string = key + "\t" + np.array2string(counts[key], separator='\t')[1:-1] + "\n"
        fh_out.write(write_string)

print("Not 20 bp:",not_20_bp)