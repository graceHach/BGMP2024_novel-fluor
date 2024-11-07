#!/usr/bin/env python

import argparse
import gzip
import glob
import numpy as np

# Run using helper script generate_counts.sh

# input directory, output directory 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--in_dir",type=str,help="directory containing merged, primer trimmed reads", default="../test_fastq/mini_trimmed_data/")
    parser.add_argument("-t","--temp_file", type=str, help="file containing all UMIs, created nu generate_counts.sh", default="all_bc_mini.txt")
    parser.add_argument("-o","--out_file",type=str,help="", default="../reports/mini_blue_counts.tsv")
    return parser.parse_args()

args = get_args()

not_20_bp = 0
# read in intermmediate file, initialize the dictionary
counts = {} 
# keys are UMIs, values are numpy arrays with nine entries
with open(args.temp_file) as fh:
    for line in fh:
        # accounting for newline character
        if len(line)==21:
            counts[line[:-1]] = np.zeros(9, dtype="int8")
        else:
            not_20_bp += 1

# will entries come out the same way they come in?
help=[]
# read in files, hash to a dictionary keyed by sequence
files = glob.glob(args.in_dir+"*.fastq.gz")
for file_index, file in enumerate(files):
    with gzip.open(file, 'rt') as fh:
        for index, line in enumerate(fh):
            if index%4==1:
                help.append(line[:-1])
                if line[:-1] in counts:
                    counts[line[:-1]][file_index] += 1

# write list of file names as header of file
header_line = "Barcode\t" + '\t'.join(files) + "\n"

keys = counts.keys()
with open(args.out_file, 'wt') as fh_out:
    fh_out.write(header_line)
    for key in keys:
        write_string = key + "\t" + np.array2string(counts[key], separator='\t')[1:-1] + "\n"
        fh_out.write(write_string)

# 