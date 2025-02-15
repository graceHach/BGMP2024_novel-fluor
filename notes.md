# 2025-02-04

Created conda environment for running nextflow with bbmerge + hts primers.

```
conda create -p countbins
mamba activate countbins/
mamba install -y nextflow 
mamba install bioconda::bbmap
mamba install bioconda::htstream
mamba deactivate
```

Run the above commands in the top of this repo. Use `mamba activate countbins/` before running `main.nf`.

I moved the nextflow script to the top directory. This will work differently than I think Grace had originally intended to run it. The pipeline will work on a single dataset at a time. 

Run one nextflow pipeline for each "red" and "blue" dataset. 

# 2025-02-15

## Pipeline

Using default parameters, for future use change parameters: 
- `params.rawfastqs_R12`
- `forward_primers`
- `reverse_primers`
 
```
(/gpfs/projects/bgmp/shared/groups/2024/novel-fluor/wesg/BGMP2024_novel-fluor/countbins) [wesg@n0349|main|BGMP2024_novel-fluor]$ nextflow main.nf

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [mad_bhabha] DSL2 - revision: 610524a1e7

executor >  local (19)
[e7/c5e2ab] merge_reads (4) | 9 of 9 ✔
[52/3a2299] trim_reads (9)  | 9 of 9 ✔
[c5/c68981] generate_counts | 1 of 1 ✔
Completed at: 14-Feb-2025 09:52:20
Duration    : 1d 19h 56m 10s
CPU hours   : 139.0
Succeeded   : 19

```
