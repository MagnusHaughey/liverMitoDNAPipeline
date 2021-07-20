

file=$1
			
file_no_extension=$(echo $file | sed 's|.*/||')	
file_no_extension=${file_no_extension::-4}

IFS='_' # underscore is set as delimiter
read -ra params <<< "$file_no_extension" # Read file parameters separately (e.g. patient, sample, primer, etc...) from file name
	
IFS=' '
unzip -c $file $file_no_extension/summary.txt | grep FAIL | cut -f 1-2
unzip -c $file $file_no_extension/summary.txt | grep WARN | cut -f 1-2

