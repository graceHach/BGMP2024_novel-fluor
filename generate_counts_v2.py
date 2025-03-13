#!/usr/bin/env python

import argparse
import gzip
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--file_list", nargs='+', type=str)
    parser.add_argument("-o","--out_file", type=str, help="", default="../reports/blue_counts_better_script.tsv")
    return parser.parse_args()

args = get_args()
files = args.file_list # files == ['file2 file3 file1 file4']
files = sorted(files[0].split()) # files == ['file1', 'file2', 'file3', 'file4']

not_20_bp = 0
print("initializing counts dict", flush=True)
counts = {} 
# keys are UMIs, values are numpy arrays with num_bins entries

for file_index, file in enumerate(files):
    print("processing: ",file, flush=True)
    with gzip.open(file, 'rt') as fh:
        for index, line in enumerate(fh):
            if index%4==1:
                if len(line)==21:
                    if line[:-1] in counts:
                        counts[line[:-1]][file_index] += 1
                    else:
                        counts[line[:-1]] = np.zeros(9, dtype="int64")
                        counts[line[:-1]][file_index] = 1
                else:
                    not_20_bp += 1

# write list of file names as header of file
header_line = "Barcode\t" + '\t'.join(files) + "\n"

print("Starting writing")
keys = counts.keys()
with open(args.out_file, 'wt') as fh_out:
    fh_out.write(header_line)
    for key in keys:
        write_string = key + "\t" + np.array2string(counts[key], separator='\t')[1:-1] + "\n"
        fh_out.write(write_string)

print("Not 20 bp:",not_20_bp)
