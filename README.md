# Per-Bin UMI Counts for Fluorescence Inference 
 
Grace Hach, Mahmoud al Mahmoud, Wesley Gomersall

Created for Dr. Calin Plesa

## Nextflow Inputs: 

Address to dir containing binned Illumina data in `fastq.gz` format.
This address must specify paired-end data as such: `*_R{1,2}*.fastq.gz`

Must also specify primer sequences for removal with htseq

## Dependencies: 

| Software | Version | 
| --- | --- | 
| bbmap | 39.15 |
| htstream | 1.4.1 | 
| nextflow | 24.10.4 | 

# User Guide

This pipeline processes paired-end Illumina data, determines the bin distribution, then integrates the PacBio data and FACS calibration 
data in order to estimate the fluorescence of novel protein sequences.
It is designed to be used after !(this repo)[https://github.com/wesleygomersall/BGMP2024_novel-fluor/edit/main/notes.md], 
which processes the corresponding PacBio data, comprising the variant and UMI sequences.

## Setup

Create the following conda environment by running these commands. Do this in the top level directory of the repo:
```
conda create -p countbins
mamba activate countbins/
mamba install -y nextflow 
mamba install bioconda::bbmap
mamba install bioconda::htstream
```

## Running main.nf

With conda environment countbins/ activated, run: 
```
nextflow main.nf --color_label "<COLOR>" --infq "<DIRECTORY CONTAINING BLUE ILLUMINA DATA>/*_R{1,2}_001.fastq.gz"
```

Note that this defaults to the primers provided along with this repo in the primers/ directory. To specify alternate primer files, run:
```
nextflow main.nf --color_label "<COLOR>" --infq "<DIRECTORY CONTAINING BLUE ILLUMINA DATA>/*_R{1,2}_001.fastq.gz" --fwd "<FORWARD PRIMER FILE>" --rev "<REVERSE PRIMER FILE>"
```
This will output the file(s) <COLOR>_counts.tsv in the output directrory.

## Running Fluorescence Estimation

The script complete_fluorescence_estimation.py integrates the two datasets and produces the final output: A csv file containing the protein 
sequences, bin distribution in all colors, and estimated fluorescence in all colors. The script is written to use an arbitrary number of bins 
and colors, and thus takes lists of arguments, enclosed in quotes, separated by spaces. For example:

```
--argument "item_1 item_2 etc"
```

It takes the following arguments:
```
-if, --illumina_counts_files (required)
Paths to one or more Illumina counts files in TSV format, as output by main.nf. The number of files provided must match the number of colors.
Format: Space-separated list of file paths.
Default: "../outputs/blue_counts.tsv ../outputs/red_counts.tsv"

-p, --prot_csv_paths (required)
Paths to input CSV files from the PacBio pipeline. Each file should contain barcode information (column 1), variant sequences (column 3), and protein sequences (column 4). The number of files must match the number of Illumina counts files and colors.
Format: Space-separated list of file paths.
Default: "../../../shared/dat/NF_results/BLUE/09_final_output/ ../../../shared/dat/NF_results/RED/09_final_output/"
(Note: Trailing slashes in the paths are necessary.)

-c, --color_list (required)
A space-separated list of colors corresponding to illumina_counts_files and prot_csv_paths. Used for labeling data.
Default: "blue red"

-o, --out_file (optional)
Path to the output CSV file where processed data will be saved.
Default: "../reports/protein_counts_and_MEFL_FINAL.csv"

-b, --bin_calibration (required)
Paths to TSV files containing bin-to-MEFL calibration data. Each file must have three columns:

    Column 1: Bin (unused)
    Column 2: Min MEFL
    Column 3: Max MEFL
    The number of files must match color_list and prot_csv_paths, and the order must be consistent.
    Format: Space-separated list of file paths.
    Default: "../bins/blue_bins.tsv ../bins/red_bins.tsv"
```

## Troubleshooting

If issues are encountered or clarification is needed, please reach out to me at ghach@uoregon.edu for support.
