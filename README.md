# NGS analysis workflow for Passman et al 2021

This repository contains all scripts necessary to implement mtDNA-seq analysis from Passman et al 2021. Input data is available via the European Genome Archive (EGA), accession number ...... 

This analysis was impemented under the Snakemake framework and executed on a HPC using the SGE system. Developed using Snakemake v6.5.3, python v3.8.5, R v3.6.1, FastQC v0.11.9, samtools v1.9. 

Prior to running the pipeline, file paths in config.yaml must be changed. Execute a "dry run" of the pipeline using the command

```
snakemake -n
``` 

To run the full pipeline, execute 

```
bash SGEclusterSubmit.sh
```

## Other useful scripts

The following scripts perform various helpful tasks on the mtDNA somatic variant calls, which are summarised in files ending in "...summary.dat" and have the format: variant\_position , frequency (of variant in test minus control) , p-value. A description of each script and how to execute them follows.

1. First sort all somatic variant summary data files into directories for each patch & patient by executing sort\_files.sh

```
bash ./sort_files.sh [path]
```

where {path} is the file path to the results/SNVs/ directory where all of the "...summary.dat" files are output by Snakemake.


2. SNV\_checkReplicates.py; to generate a plot of variant frequencies in replicate 1 vs. replicate 2

```
python3 ./SNV_checkReplicates.py [--replicateA A] [--replicateB B] [--output O]
```

where\
&nbsp;  --replicateA A &emsp;&emsp;	path to file with A replicate data (ending in "...repA...summary.dat")\
&nbsp;  --replicateB B &emsp;&emsp;     path to file with corresponding B replicate data (ending in "...repB...summary.dat")\
&nbsp;  --output O &emsp;&emsp;     path to output figure (will be appended with "...\_replicate\_frequencies.png")\


3. To run the above script on a directory of files ending in "...\_summary.dat", execute the compare\_replicates.sh shell script.

```
bash ./compare_replicates.sh [path]
```

where {path} is the file path to the directory containing "...summary.dat" files.


4. compile\_all\_plots.py; to produce a .tex file with a page for each patient+sample containing output deepSNV scatter plot for both replicate, and the replicate comparison figure produced with the SNV\_checkReplicates.py script. Requires directory containing the deepSNV scatter plots for each patient+sample (files ending in "...scatterPlot.pdf"). This should be the same directory as that containing the "...\_replicate\_frequencies.png" files.

```
python3 ./compile_all_plots.py [path]
```

where {path} is the file path to directory containing "...scatterPlot.pdf" and "...\_replicate\_frequencies.png" files.


5. shared\_and\_private\_variants.py; sorts shared (public) and private variants for each patient and summarises these in a .dat file.

```
python3 shared_and_private_variants.py [path]
```

where {path} is the file path to the directory containing "...summary.dat" files. Must have sorted somatic calls using {./sort\_files.sh} prior to executing this script.


