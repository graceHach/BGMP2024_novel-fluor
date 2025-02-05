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

