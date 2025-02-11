#!/usr/bin/env python
import argparse
import glob

def get_args():
    parser= argparse.ArgumentParser()
    # parser.add_argument('-cs','--column_sums', nargs='+', help='Columns sums from the tsv', required=True)
    parser.add_argument('-f','--in_tsv', type=str, help='Input tsv file, barcodes (first column) and counts in all bins. Has header', 
                        default="../../../shared/dat/counts_per_bin/blue_counts_FINAL.tsv")
    parser.add_argument("-o",'--out_file', type=str, help='output tsv file',default="../../../shared/dat/counts_per_bin/blue_prot_counts.tsv")
    parser.add_argument("-p",'--prot_csv_path', type=str, 
                        help='Path to input csv files, 1st column barcodes, 3rd column variant sequences, 4th column protein sequences',
                        default="../../../shared/dat/NF_pacbio_output/blu/09_final_output/")
    return parser.parse_args()

args = get_args()
num_bins = 9 # len(args.column_sums)
# col_sums = [int(x) for x in args.column_sums]
# grand_total = sum(col_sums)
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

first = True
with open(args.out_file, 'wt') as fo, open(args.in_tsv, 'rt') as fin:
    for line in fin:
        if first:
            first=False
            header=line # header contains file name
            new_header="protein\t"+"\t".join(["bin"+str(i) for i in range(1,num_bins+1)])+"\ttotalreads\n"
            fo.write(new_header)
        else:
            split = line.split(sep="\t")
            counts_list = [int(x) for x in split[1:num_bins+1]]
            # norm_counts_list = list(map(lambda x,y: str(x/y), counts_list,col_sums))
            total_number_of_reads = sum(counts_list)
            # noramalized_total = total_number_of_reads/grand_total
            barcode = split[0]
            if barcode in seq_dict:
                # write_line = seq_dict[barcode]+"\t"+"\t".join(norm_counts_list)+"\t"+str(total_number_of_reads)+"\t"+str(noramalized_total)+"\n"
                write_line = seq_dict[barcode]+"\t"+"\t".join(map(str, counts_list))+"\t"+str(total_number_of_reads)+"\n"
                fo.write(write_line)
