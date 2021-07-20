# NGS analysis workflow for Passman et al 2021

This repository contains all scripts necessary to implement mtDNA-seq analysis from Passman et al 2021. Input data is available via the European Genoma Archive (EGA), accession number ...... 

This analysis was impemented under the Snakemake framework and executed on a HPC using the SGE system. Developed using Snakemake v6.5.3, python v3.8.5, R v3.6.1, FastQC v0.11.9, samtools v1.9. 

Prior to running the pipeline, file paths in config.yaml must be changed. Execute a "dry run" of the pipeline using the command

```
snakemake -n
``` 

To run the full pipeline, execute 

```
bash SGEclusterSubmit.sh
```

## Other useful scripts
