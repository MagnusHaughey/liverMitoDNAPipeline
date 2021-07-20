

#/bin/bash 

for path in $1/*A_summary.dat;
do

	file=$(echo $path | rev | cut -d "/" -f 1 | rev)


	IFS='_'
	read -ra fragments <<< "$file"

	IFS=' '
	primer=$(echo ${fragments[5]} | cut -c1-2)

	filepath=
	basic_filename=$1'/'${fragments[0]}'_'${fragments[1]}'_'${fragments[2]}'_'${fragments[3]}'_'${fragments[4]}'_'$primer


	#=========================	For germline removed data
	replicateA=$basic_filename'A_summary.dat'
	replicateB=$basic_filename'B_summary.dat'
	outfile_name=$basic_filename'_shared_variants.dat'
	


	if [[ $basic_filename != *"Stroma"* ]] && [[ ! -f $outfile_name ]] && [ -f "$replicateA" ] && [ -f "$replicateB" ]; then

		printf "Processing '$basic_filename'"

		# Use absolute path to "SNV_checkReplicates.py"
		#python3 ./absolute_path_to/SNV_checkReplicates.py \
		#python3 /data/BCI-EvoCa2/magnus/liver_mtDNA/mtDNApipeline/liverMitoDNAPipeline/otherUsefulScripts/SNV_checkReplicates.py \
		python3 ./SNV_checkReplicates.py \
				--replicateA $replicateA \
				--replicateB $replicateB \
				--output $basic_filename

	fi
	printf "\n"

done
