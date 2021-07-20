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

The following scripts perform various helpful tasks on the mtDNA somatic variant calls, which are summarised in files ending in "...\_summary.dat" and have the format: variant\_position , frequency (of variant in test minus control) , p-value. A description of each script and how to execute them follows.

1) First sort all variant summary data files into directories for each patch & patient by executing sort\_files.sh

	Usage:
		bash sort_files.sh [arguments]
	
	Arguments:
		Path to directory containing all SNV data output from Snakemake pipeline


2) SNV\_checkReplicates.py; to generate a plot of variant frequencies in replicate 1 vs. replicate 2

	Usage:
		python3 SNV_checkReplicates.py [arguments]
	
	Arguments:
		--replicateA	path to file with A replicate data (ending in "..._summary.dat")
		--replicateB	path to file with B replicate data (ending in "..._summary.dat")
		--output	path to output figure (will be appended with "..._replicate_frequencies.png")


3) To run the above script on a directory of files ending in "...\_summary.dat", execute the compare\_replicates.sh shell script.

	Usage:
		bash compare_replicates.sh [arguments]
	
	Arguments:
		Path to directory containing "...summary.dat" files.


4) compile\_all\_plots.py; to produce a .tex file with a page for each patient+sample containing output deepSNV scatter plot
   for both replicate, and the replicate comparison figure produced with the SNV\_checkReplicates.py script. Requires directory
   containing the deepSNV scatter plots for each patient+sample (files ending in "...scatterPlot.pdf"). This should be the
   same directory as that containing the "...\_replicate\_frequencies.png" files.

	Usage:
		python3 compile_all_plots.py [arguments]
	
	Arguments:
		 Path to directory containing "...scatterPlot.pdf" and "..._replicate_frequencies.png" files.

5) shared\_and\_private\_variants.py; sorts shared and private variants for each patient and summarises these in a .dat file.

	Usage:
		python3 shared_and_private_variants.py [arguments]

	Arguments:
		path to file containing "..._summary.dat" data files. File path should specify the parent directory AND the generic 
		beginning of the file name e.g. "/path_to_parent_directory/GC_EC_8466_193B"

	Output:
		`{prefix}_rep1.dat`
		`{prefix}_rep2.dat`

