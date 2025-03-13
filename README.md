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
It is designed to be used after (this repo)[https://github.com/wesleygomersall/BGMP2024_novel-fluor/edit/main/notes.md], 
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
and colors. It takes the following arguments:

```

```
