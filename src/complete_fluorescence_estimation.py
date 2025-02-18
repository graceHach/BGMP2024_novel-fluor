#!/usr/bin/env python

import argparse
import gzip
import glob
import numpy as np
import pandas as pd

'''
Downstream of generate_counts_v2.py
Replaces generate_protein_counts.py, counts_matrix_manipulation.py, and bins_to_MEF.py.
Replaces estimate_fluorescence.sh, with fully python implementation
'''

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-if","--illumina_counts_files", nargs='+', type=str, default=["../reports/blue_counts.tsv ../reports/red_counts.tsv"],
                        help="Paths to Illumina counts files in tsv format as output by generate_counts_v2.py. Must be as many counts files as \
                            there are colors. Format as space-separated list (see default)")
    parser.add_argument("-p",'--prot_csv_paths', nargs='+', type=str, 
                        help='Paths to input csv files as space separated list. 1st column barcodes, 3rd column variant sequences, 4th column \
                        protein sequence. must be as many paths as there are count files/colors above.',
                        default=["../../../shared/dat/NF_results/BLUE/09_final_output/ ../../../shared/dat/NF_results/RED/09_final_output/"])
                        # trailing slashes are totally necessary in above paths, IDK why
    parser.add_argument("-c",'--color_list', nargs='+', type=str, 
                        help='space-separated list of colors corresponding to illumina_counts_files/prot_csv_paths. Used for labeling dataframes',
                        default=["blue red"])
    parser.add_argument("-o","--out_file", type=str, help="", default="../reports/complete_fluorescence_estimation.csv")
    return parser.parse_args()

args = get_args()


############################### CONVERT BARCODES FROM COUNTS FILES TO PROTEINS ############################### 

# read in illumina counts file paths
illu_files = args.illumina_counts_files # formatted as ['file2 file3 file1 file4']
# In same order as prot_csv_paths, color_list
illu_files = illu_files[0].split() # formatted as ['file1', 'file2', 'file3', 'file4']
color_list = args.color_list
color_list = color_list[0].split()
num_colors = len(color_list)


# Reading first count file
df = pd.read_csv(illu_files[0], sep="\t")

# Change column names to be "barcode, <color> 1, <color> 1, <color> 3, etc"
num_bins = len(df.columns)-1
new_col_names = ['barcode'] + [color_list[0]+" "+str(x+1) for x in range(num_bins)]
color_list = color_list[1:] # remove the item from color list after it's been used
df.columns = new_col_names

# Reading in subsequent  counts files
illu_files = illu_files[1:]

for index, file in enumerate(illu_files):
    new_df = pd.read_csv(file, sep="\t")
    num_bins = len(new_df.columns)-1
    new_col_names = ['barcode'] + [color_list[index]+" "+str(x+1) for x in range(num_bins)]
    new_df.columns = new_col_names
    df = df.merge(new_df, how="outer", on='barcode')

# Now we have a dataframe with the following columns (where n is the number of bins and m is the number of colors)
# barcode, <color 1> 1, ..., <color 1> n_bins, <color 2> 1, ..., <color 2> n_bins, ..., ..., <color m> n_bins

############# Read in protein data
prot_paths = args.prot_csv_paths # formatted as ['file2 file3 file1 file4']
# In same order as prot_csv_paths, color_list
prot_paths = prot_paths[0].split() # formatted as ['file1', 'file2', 'file3', 'file4']

# Get a dict of all barcodes from pandas dataframe
bc_to_prot = {key:"" for key in df['barcode']}
# Initialize each barcode to an empty string

# Conflics between colors and bins must be resolved
# How?
conflicts = 0
not_in_there = 0
# Note that each path is a color
for path in prot_paths:
    files = sorted(list(glob.glob(path+"*")))
    # Each file is a bin. They are zero-indexed (in terms of bin labels) and in csv format.
    for file in files:
        # Column 0 is bc, column 3 is protein sequence
        prot_df = pd.read_csv(file)
        # Populate dict with protein sequences
        barcodes, proteins = list(prot_df.iloc[:,0]), list(prot_df.iloc[:,3])
        for index, bc in enumerate(barcodes):
            if bc in bc_to_prot:
                # If in there and empty, overwrite
                if bc_to_prot[bc] == "":
                    bc_to_prot[bc] = proteins[index]
                # If the barcode is in the dict and already filled, don't overwrite.
                # Earlier bins are larger, so barcode-bin assignments on lower integer
                # bins are presumed to be the correct ones.
                elif not bc_to_prot[bc] == proteins[index]:
                    # If the new sequence assignment conflics with the old, this is a conflict
                    conflicts +=1
            else:
                # count the PacBio barcodes that are not in the illumina dictionary
                # I'm not adding these to the dict because there are no counts, so fluorescence estimation 
                # is meaningless
                not_in_there += 1


num_empty_prots = sum(1 for value in bc_to_prot.values() if value == '')
print("conflics:", conflicts)
print("not_in_there:", not_in_there)
print("populated:", len(bc_to_prot))
print("still empty:", num_empty_prots)
