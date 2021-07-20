

To process deepSNV variant data which is outputted from the mtDNA pipeline, execute the following scripts. In the following, 
files with filename ending in "..._summary.dat" contain re-formatted deepSNV variant data, and the data within are arranged in
the following format: variant_position , frequency (of variant in test minus control) , p-value.

1) First sort all variant summary data files into directories for each patch & patient by executing sort_files.sh

	Usage:
		bash sort_files.sh [arguments]
	
	Arguments:
		Path to directory containing all SNV data output from Snakemake pipeline


2) SNV_checkReplicates.py; to generate a plot of variant frequencies in replicate 1 vs. replicate 2

	Usage:
		python3 SNV_checkReplicates.py [arguments]
	
	Arguments:
		--replicateA	path to file with A replicate data (ending in "..._summary.dat")
		--replicateB	path to file with B replicate data (ending in "..._summary.dat")
		--output	path to output figure (will be appended with "..._replicate_frequencies.png")


3) To run the above script on a directory of files ending in "..._summary.dat", execute the compare_replicates.sh shell script.

	Usage:
		bash compare_replicates.sh [arguments]
	
	Arguments:
		Path to directory containing "...summary.dat.noGermline.dat" files.


4) compile_all_plots.py; to produce a .tex file with a page for each patient+sample containing output deepSNV scatter plot
   for both replicate, and the replicate comparison figure produced with the SNV_checkReplicates.py script. Requires directory
   containing the deepSNV scatter plots for each patient+sample (files ending in "...scatterPlot.pdf"). This should be the
   *SAME DIRECTORY* as that containing the "..._replicate_frequencies.png" files.

	Usage:
		python3 compile_all_plots.py [arguments]
	
	Arguments:
		 Path to directory containing "...scatterPlot.pdf" and "..._replicate_frequencies.png" files.

5) shared_and_private_variants.py; sorts shared and private variants for each patient and summarises these in a .dat file.

	Usage:
		python3 shared_and_private_variants.py [arguments]

	Arguments:
		path to file containing "..._summary.dat" data files. File path should specify the parent directory AND the generic 
		beginning of the file name e.g. "/path_to_parent_directory/GC_EC_8466_193B"

	Output:
		`{prefix}_rep1.dat`
		`{prefix}_rep2.dat`









