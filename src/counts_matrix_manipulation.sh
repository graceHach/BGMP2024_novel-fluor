#!/bin/bash

echo "The truth will set you free. But not until it is finished with you."

in_file="../../../shared/dat/counts_per_bin/blue_prot_counts.tsv"
out_file="../../../shared/dat/counts_per_bin/blue_median_bin.tsv"

./counts_matrix_manipulation.py --in_tsv $in_file --out_file $out_file

echo "finished blue"

in_file="../../../shared/dat/counts_per_bin/red_prot_counts.tsv"
out_file="../../../shared/dat/counts_per_bin/red_median_bin.tsv"

./counts_matrix_manipulation.py --in_tsv $in_file --out_file $out_file

echo "finished red"