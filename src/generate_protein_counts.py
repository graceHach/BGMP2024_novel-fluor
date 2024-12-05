#!/usr/bin/env python3
import argparse
import glob

def get_args():
    parser= argparse.ArgumentParser()
    parser.add_argument('-cs','--column_sums', nargs='+', help='Columns sums from the tsv', required=True)
    parser.add_argument('-f','--in_tsv', type=str, help='Input tsv file, barcodes (first column) and counts in all bins. Has header', 
                        default="../../../shared/dat/counts_per_bin/blue_counts_FINAL.tsv")
    parser.add_argument("-o",'--out_file', type=str, help='output tsv file',default="../../../shared/dat/counts_per_bin/blue_prot_counts.tsv")
    parser.add_argument("-p",'--prot_csv_path', type=str, 
                        help='Path to input csv files, 1st column barcodes, 3rd column variant sequences, 4th column protein sequences',
                        default="../../../shared/dat/NF_pacbio_output/blu/09_final_output/")
    return parser.parse_args()

args = get_args()
num_bins = len(args.column_sums)
col_sums = [int(x) for x in args.column_sums]
grand_total = sum(col_sums)
# Populate dictionary of barcodes (keys) and protein sequences (values) from each file in prot_csv_path
seq_dict = {}
files = sorted(list(glob.glob(args.prot_csv_path+"*")))
for file in files:
    print(file)
    with open(file, 'rt') as fh:
        for line in fh:
            line=line.split(sep=",")
            # column 0 is bc
            # column 2 is dna seq
            # column 3 is protein sequence
            if line[0] not in seq_dict:
                if not line[3] == '':
                    seq_dict[line[0]] = line[3]

# Re-written to store each line in a dictionary, keyed to the aa sequence. If a duplciate aa sequence is reached,
# the rows are added. This is fine because the normalized count list and the normalized total are fractions with a 
# common denominator.

prot_dict = {}
# key is aa sequence, value is a list with the elements:
# [normalized_counts_list, total_number_of_reads, normalized_total]
first = True
with open(args.in_tsv, 'rt') as fin:
    for line in fin:
        if first:
            first=False
        else:
            split = line.split(sep="\t")
            # Get a list of counts for that UMI
            counts_list = [int(x) for x in split[1:num_bins+1]]
            # Normalize the counts list by the total sums of each column
            # Keep these as floats and not strings because they may need to be summed later
            norm_counts_list = list(map(lambda x,y: x/y, counts_list,col_sums))
            # Get the total number of each row
            total_number_of_reads = sum(counts_list)
            noramalized_total = total_number_of_reads/grand_total
            barcode = split[0]
            if barcode in seq_dict:
                protein_sequence = seq_dict[barcode]
                if protein_sequence in prot_dict:
                    # Duplicate encountered. Combine all entries by summing:
                    # Combine normalized_counts_list
                    prot_dict[protein_sequence][0] = [x+y for x,y in zip(prot_dict[protein_sequence][0], norm_counts_list)]
                    # Combine total number of reads
                    prot_dict[protein_sequence][1] = prot_dict[protein_sequence][1] + total_number_of_reads
                    # Combine normalized total
                    prot_dict[protein_sequence][2] = prot_dict[protein_sequence][2] + noramalized_total 
                else:
                    # Not a duplicate sequence (so far)
                    prot_dict[protein_sequence] = [norm_counts_list, total_number_of_reads, noramalized_total]

# Use dict to write out to file
first = True
with open(args.out_file, 'wt') as fo:
    if first:
        new_header="protein\t"+"\t".join(["bin"+str(i) for i in range(1,num_bins+1)])+"\ttotalreads\ttotalfrac\n"
        fo.write(new_header)
    for protein_sequence in prot_dict:
        norm_counts_list, total_number_of_reads, noramalized_total = prot_dict[protein_sequence]
        # Needs to be converted back to string
        norm_counts_list = [str(x) for x in norm_counts_list]
        write_line = protein_sequence+"\t"+"\t".join(norm_counts_list)+"\t"+str(total_number_of_reads)+"\t"+str(noramalized_total)+"\n"
        fo.write(write_line)
