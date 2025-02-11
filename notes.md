2025-02-04

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
