#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def get_args():
    parser= argparse.ArgumentParser()
    parser.add_argument('-f','--in_tsv', type=str, help='Input tsv file, barcodes (first column) and counts in all bins. Has header', 
                        default="../../../shared/dat/counts_per_bin/blue_prot_counts.tsv")
    parser.add_argument("-o",'--out_file', type=str, help='output tsv file',default="../../../shared/dat/counts_per_bin/blue_matrix.tsv")
    return parser.parse_args()

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

args = get_args()

# read in dataframe
# columns: 0 protein	1 bin1	2 bin2	3 bin3	4 bin4	5 bin5	6 bin6	7 bin7	8 bin8	9 bin9	10 totalreads	11 totalfrac
in_df = pd.read_csv(args.in_tsv, sep='\t')
n_bins = in_df.shape[1] - 3
counts_columns = list(in_df.columns)[1:n_bins+1]
# add fold change columns by normalizing using totalfrac
for col in counts_columns:
    in_df[col+"_fc"] = in_df[col]/in_df["totalfrac"]
print(in_df.columns)

# add A matrix by normalizing each fc column by rowsum of all fc columns
fc_columns = sorted(list(in_df.columns)[n_bins+3:])
print(fc_columns)
# axis=1 sums across rows 
in_df["rowsum"] = in_df[fc_columns].sum(axis=1)
for index, fc_col in enumerate(fc_columns):
    in_df["A"+str(index+1)] = in_df[fc_col]/in_df['rowsum']
print(in_df.columns)

# Create new dataframe proteinAcs, or the A columns, but as cumulative sums of current and prev columns
cumulative_As = pd.DataFrame()
A_cols = ["A"+str(i) for i in range(1,n_bins+1)]
for i in range(1,n_bins+1):
    print(A_cols[0:i])
    cumulative_As["A_"+str(i)+"_cs"] = in_df[A_cols[0:i]].sum(axis=1)

# get vectors of the row indices that bookend 0.5 
#  not a typo
lower_indices = maxcs(cumulative_As)
upper_indices = mincs(cumulative_As)

# modify in_df to add lower_index, upper_index, median column
# is cumulative p the same as cumulative_As?
# lower_indices and upper_indices have None values 
print(cumulative_As)
print(sum([1 for x in lower_indices if not x]))
print(len(lower_indices))
#print(sum(np.isnan(upper_indices)))

in_df["lower_index"] = lower_indices
in_df["upper_index"] = upper_indices

# ONLY THE MEDIAN COLUMN IS USED IN THE FINAL CALCULATION
out_df = pd.DataFrame()
out_df["protein"] = in_df["protein"]
out_df["median"]= calc_median(cumulative_As,upper_indices,lower_indices)

out_df.to_csv(args.out_file, index=False)

print(out_df)