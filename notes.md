2025-02-04

Created conda environment for running nextflow with bbmerge + hts primers.

```
conda create -p countbins
mamba activate countbins/
mamba install -y nextflow 
mamba install bioconda::bbmap
mamba install bioconda::htseq
mamba deactivate
```

Run the above commands in the top of this repo. Use `mamba activate countbins/` before running `main.nf`.
