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

# Running this pipeline

The following commands do not specify the `--fwd` and `--rev` primer files. Either edit the existing files to contain accurate primer sequences or specify with the fwd and rev options. 

```
$ nextflow main.nf --infq "../../shared/rawdata/BLUE/NovaSeq_GC3F_7125/blue-sort_bin*_R{1,2}_001.fastq.gz" --ou
t blue_out/

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [happy_pauling] DSL2 - revision: 80689b1676

[ac/a10bdb] merge_reads (3) | 9 of 9 ✔
[5a/85de6c] trim_reads (8)  | 9 of 9 ✔
[04/592688] generate_counts | 1 of 1 ✔
Completed at: 16-Feb-2025 21:51:38
Duration    : 1d 10h 40m 54s
CPU hours   : 197.9
Succeeded   : 19
```

and

```
$ nextflow main.nf --infq '../../shared/rawdata/RED/NovaSeq_GC3F_7124/red-sort_bin*_R{1,2}_001.fastq.gz' --out red_out/

 N E X T F L O W   ~  version 24.10.4

Launching `main.nf` [exotic_visvesvaraya] DSL2 - revision: 80689b1676

executor >  local (9)
[1b/e0b2f9] merge_reads (2) | 2 of 9
executor >  local (11)
[1b/e0b2f9] merge_reads (2) | 2 of 9
executor >  local (11)
[1b/e0b2f9] merge_reads (2) | 2 of 9
executor >  local (11)
[1b/e0b2f9] merge_reads (2) | 2 of 9
[99/f9a077] merge_reads (7) | 3 of 9
[a6/9181a6] trim_reads (3)  | 2 of 3
[99/f9a077] merge_reads (7) | 3 of 9
[a6/9181a6] trim_reads (3)  | 2 of 3
[63/50f7f0] merge_reads (4) | 9 of 9 ✔
[3e/291a75] trim_reads (9)  | 9 of 9 ✔
[03/40d99b] generate_counts | 1 of 1 ✔
Completed at: 17-Feb-2025 06:56:24
Duration    : 1d 19h 49m 23s
CPU hours   : 139.0
Succeeded   : 19
```

