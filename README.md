# Per-Bin UMI Counts for Fluorescence Inference 
 
Grace Hach, Mahmoud al Mahmoud, Wesley Gomersall

Created for Dr. Calin Plesa

## Nextflow Inputs: 

`--infq`: String enclosed in quotes (single or double) containing address to dir containing binned Illumina data in `fastq.gz` format.
This address must specify paired-end data as such: `'dir/*_R{1,2}*.fastq.gz'`

`--fwd` and `--rev`: Specify forward and revers primer sequences respectively for removal with htseq.

`--out`: output directory to publish merged and trimmed outputs, and result counts file `counts.tsv`.

## Dependencies: 

| Software | Version | 
| --- | --- | 
| bbmap | 39.15 |
| htstream | 1.4.1 | 
| nextflow | 24.10.4 | 

