#!/bin/bash

# Takes in file with bin edges and another with median bins for each variant and converts to 

# Run blue
./bins_to_MEF.py --in_bin_edges "../bins/blue_bins.tsv"  --in_median_bins "../../../shared/dat/counts_per_bin/blue_median_bin.tsv" --out_file "../reports/blue_FINAL.tsv"

# Run red
./bins_to_MEF.py --in_bin_edges  "../bins/red_bins.tsv" --in_median_bins "../../../shared/dat/counts_per_bin/red_median_bin.tsv" --out_file "../reports/red_FINAL.tsv"
