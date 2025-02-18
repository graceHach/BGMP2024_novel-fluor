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
df.columns = new_col_names

# Reading in subsequent  counts files
illu_files = illu_files[1:]

for index, file in enumerate(illu_files):
    new_df = pd.read_csv(file, sep="\t")
    num_bins = len(new_df.columns)-1
    # color list index + 1 because index of illu_files starts at zero, but color_list starts at one
    new_col_names = ['barcode'] + [color_list[index+1]+" "+str(x+1) for x in range(num_bins)]
    new_df.columns = new_col_names
    df = df.merge(new_df, how="outer", on='barcode')

# Now we have a dataframe with the following columns (where n is the number of bins and m is the number of colors)
# barcode, <color 1> 1, ..., <color 1> n_bins, <color 2> 1, ..., <color 2> n_bins, ..., ..., <color m> n_bins

############# Read in protein data
prot_paths = args.prot_csv_paths # formatted as ['file2 file3 file1 file4']
# In same order as prot_csv_paths, color_list
prot_paths = prot_paths[0].split() # formatted as ['file1', 'file2', 'file3', 'file4']

# Get a dict of all barcodes from pandas dataframe
#blu_bc_to_prot = {key:"" for key in df['barcode']}
illumina_bcs = set(df['barcode'])
# Initialize each barcode to an empty string

# Conflics between colors and bins must be resolved
# How?
between_bin_conflicts = 0
between_color_conflicts = 0
in_pacbio_not_illu = 0
what = 0
# list of list of dicts:
# outer list is each color, inner list is each bin
list_of_dicts = []
# Note that each path is a color, each file is a bin
for color_index, path in enumerate(prot_paths):
    # add a new empty list corresponding to each color
    print("color ind", color_index)
    list_of_dicts.append([])
    files = sorted(list(glob.glob(path+"*")))
    # Each file is a bin. They are zero-indexed (in terms of bin labels) and in csv format.
    for file in files:
        print(file)
        # create new empty dict for the bin
        bc_to_prot = {}
        # Column 0 is bc, column 3 is protein sequence
        prot_df = pd.read_csv(file)
        # Populate dict with protein sequences
        barcodes, proteins = list(prot_df.iloc[:,0]), list(prot_df.iloc[:,3])
        for bc,prot in zip(barcodes, proteins):
            # If the barcode has been counted in the illumina data and the protein isn't empty
            if bc in illumina_bcs and not prot=="":
                # add to dict
                bc_to_prot[bc] = prot
            else:
                # count the PacBio barcodes that are not in the illumina barcode set
                # I'm not adding these to the dict because there are no counts, so fluorescence estimation 
                # is meaningless
                in_pacbio_not_illu += 1
        list_of_dicts[color_index].append(bc_to_prot)

# unpack list of dicts, first merging each bin 
# This will be a list of <number of colors> list of dicts
second_list_of_dicts = []
for color in list_of_dicts:
    color_dict = {}
    for bin in color:
        for bc in bin.keys():
            if bc in color_dict:
                # if the barcode is already in the color dict, it's been seen before in previous bins
                if not color_dict[bc] == bin[bc]:
                    # if the vale isn't the same, it's a conflict
                    between_bin_conflicts += 1
            else:
                color_dict[bc] = bin[bc]
    second_list_of_dicts.append(color_dict)

# Now merging the dicts for each color
final_bc_prot_dict = {}
for channel_dict in second_list_of_dicts:
    for bc in channel_dict.keys():
        if bc in final_bc_prot_dict:
            # if the barcode is already in the protein dict, it's been seen before in previous colors
            if not final_bc_prot_dict[bc] == channel_dict[bc]:
                between_color_conflicts +=1
        else:
            final_bc_prot_dict[bc] = channel_dict[bc]


print("between bin conflicts:",between_bin_conflicts)
print("between color conflicts:",between_color_conflicts)
print("in Pacbio not Illumina",in_pacbio_not_illu)
print("final number of bc-prot pairs:",len(final_bc_prot_dict))

'''
between bin conflicts: 29517
between color conflicts: 8066
in Pacbio not Illumina 90581
final number of bc-prot pairs: 33851
'''

# Create dataframe from final protein dictionary
final_df = pd.DataFrame()
# These do come out in corresponding order
final_df['barcode'] = final_bc_prot_dict.keys()
final_df['protein'] = final_bc_prot_dict.values()

# Merge the two dataframes, keeping only the barcodes in final_df
final_df = final_df.merge(df, on="barcode", how='left')
# fill NaNs with zeros for those that don't have counts in a given color
final_df.fillna(0, inplace=True)
# Columns of final_df are now:
#"barcode", "protein", "<color1> 1"
# counts of each bin in each color each 


############################### CREATE DATAFRAMES WITH NORMALIZED COUNTS ############################### 

# Note that this is doing much of what get_prot_counts did, excluding converting bc to protein sequences
# which was done above

# iterate over colors, getting norm counts columns and adding tot num reads, normalized total column
list_of_color_df = []
for color_label in color_list:
    # get column labels for a given color
    col_labels = [color_label+" "+str(x+1) for x in range(num_bins)]
    # create new dataframe for a given color
    df = pd.DataFrame()
    # get column sums for a given color
    col_sums = final_df[col_labels].sum(axis=0)
    # Get grand total for the color
    grand_total = sum(col_sums)
    # Create normalized counts list, normalizing the list of integer counts to the sums of each column
    norm_counts_label = ["norm_counts_"+str(x+1) for x in range(len(col_labels))]
    df[norm_counts_label] = final_df[col_labels]/col_sums
    # sum the total number of reads for each sequence, in a given color,  add to new df
    df["totalreads"] = final_df[col_labels].sum(axis=1)
    # Normalize total number of reads to the grand total 
    df["totalfrac"] = df["totalreads"]/grand_total
    # add to list of dfs
    list_of_color_df.append(df)

############################### MATRIX MANIPULATION ON NORMALIZED COUNTS ############################### 

# This script replaces counts_matrix_manipulation.py
# Based on preliminary r script by Andrew Holston, who is the best.

# Define functions
def maxcs(df):
    """
    takes a dataframe with only cumulative A values, returns the maximum row index where the value is less than 0.5
    Input:

    Output:
        numpy array last_indices where each entry is a row
        If there are no values < 0.5 in the cumulative dataset the dataframe contains none at the index
    
    Original function:
    maxcs = function(x, output){
    return(max(which(c(x[1],x[2],x[3],x[4],x[5],x[6]) < 0.5)))
    }
    Some of these
    """
    logical=df<0.5
    logical_arr = logical.to_numpy()
    # Initialize an array with None
    last_indices = np.full(logical_arr.shape[0], None)
    # Iterate over rows and get the largest index of "True" for each row
    for index, row in enumerate(logical_arr):
        true_indices = np.where(row)[0]  # Get indices where row is True
        if true_indices.size > 0:
            last_indices[index] = true_indices[-1]  # Get the last index
    return last_indices

def calc_median(cumulative_As, upper_index, lower_index):
    """
    Input:
        cumulative_As - numpy matrix of cumulative A values, each row is analagous to cumulative_p

    Output:

    Replaces:
    median = ifelse(
      is.infinite(lower_index), 
      1, 
      lower_index + (0.5 - unlist(cumulative_p)[lower_index]) / 
        (unlist(cumulative_p)[upper_index] - unlist(cumulative_p)[lower_index])
    )
    """
    cumulative_As = cumulative_As.to_numpy()
    medians = np.zeros(cumulative_As.shape[0])
    for index, row in enumerate(cumulative_As):
        if not lower_index[index]:
            medians[index] = 1
        else:
            medians[index] = lower_index[index] + (0.5 - row[lower_index[index]])/(row[upper_index[index]] - row[lower_index[index]])
    return medians 

def mincs(df):
    """
    takes a dataframe with only cumulative A values, returns the maximum row index where the value is greater than or equal to 0.5
    Original function:
    mincs = function(x, output){
    return(min(which(c(x[1],x[2],x[3],x[4],x[5],x[6]) >= 0.5)))
    }
    If a vector is all
    """
    logical=df>=0.5
    logical_arr = logical.to_numpy()
    # Initialize an array with None
    last_indices = np.full(logical_arr.shape[0], None)
    # Iterate over rows and get the largest index of "True" for each row
    for index, row in enumerate(logical_arr):
        true_indices = np.where(row)[0]  # Get indices where row is True
        if true_indices.size > 0:
            last_indices[index] = true_indices[0]  # Get the first index
    return last_indices

print("got this far")
# Iterate pandas dataframes, representing each color
for df in list_of_color_df:
    # Columns are:
    #       0 - norm_count_1
    #       num_bins-1 - norm_count_n_bins 
    #       num_bins - totalreads 
    #       num_bins + 1 - 'totalfrac
    # Almost the same as in counts_matrix_manipulation, but without a protein column in index 0
    #print(df.columns)
    counts_columns = list(df.columns)[:num_bins]
    # add fold change columns by normalizing using totalfrac
    for col in counts_columns:
        df[col+"_fc"] = df[col]/df["totalfrac"]

    # add A matrix by normalizing each fc column by rowsum of all fc columns
    fc_col_labels = [col +"_fc" for col in counts_columns]
    # axis=1 sums across rows 
    df["rowsum"] = df[fc_col_labels].sum(axis=1)
    for index, fc_col in enumerate(fc_col_labels):
        df["A"+str(index+1)] = df[fc_col_labels]/df['rowsum']

    '''
    # Create new dataframe proteinAcs, or the A columns, but as cumulative sums of current and prev columns
    cumulative_As = pd.DataFrame()
    A_cols = ["A"+str(i) for i in range(1,n_bins+1)]
    for i in range(1,n_bins+1):
        #print(A_cols[0:i])
        cumulative_As["A_"+str(i)+"_cs"] = in_df[A_cols[0:i]].sum(axis=1)
    '''
    

