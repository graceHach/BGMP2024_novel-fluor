#!/usr/bin/env python

import argparse
import pandas as pd

# input directory, output directory 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ie","--in_bin_edges",type=str,help="tsv file containing bin edges in MEFL with column headings: bin, bin minimum, bin maximum", default="../bins/blue_bins.tsv")
    parser.add_argument("-im","--in_median_bins",type=str,help="tsv file containing bin edges in MEFL with column headings: bin, bin minimum, bin maximum", default="../../../shared/dat/counts_per_bin/blue_median_bin.tsv")
    parser.add_argument("-o","--out_file",type=str,help="", default="../reports/blue_MEF.tsv")
    return parser.parse_args()

args = get_args()

# Header: Bin	min_MEFL	max_MEFL
in_bin_boundaries = pd.read_csv(args.in_bin_edges, sep='\t')
# Convert to list where index + 1 = bin

# Header: protein	median
in_median_bins = pd.read_csv(args.in_median_bins, sep='\t')

# Create new column with MEF
# Calcualted for a given bin as bin_minimum + (bin max - bin min) * fractional part of median bin
# Having trouble with strict pandas implementation 
min_bins_list = list(in_bin_boundaries["min_MEFL"])
max_bins_list = list(in_bin_boundaries["max_MEFL"])
medians_list = list(in_median_bins["median"])
medians_MEF = []
for median in medians_list:
    median_index = int(median) - 1
    median_MEF = min_bins_list[median_index] + (max_bins_list[median_index] - min_bins_list[median_index])*(median%1)
    medians_MEF.append(median_MEF)

#out_df = pd.DataFrame(in_median_bins["protein"], medians_MEF)
in_median_bins["MEF"] = medians_MEF
#out_df.columns = ["protein", "MEF"]
# Delete median bin column
in_median_bins.drop('median', axis=1, inplace=True)
print(in_median_bins)
in_median_bins.to_csv(args.out_file, index=False)




