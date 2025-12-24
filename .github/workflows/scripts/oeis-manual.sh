#!/bin/bash

# Ensure the factorization executable is available
if [[ ! -x "./factor" ]]; then
	echo "Error : The 'factor' executable isn't available."
	echo "Please place the executable in the current directory."
	exit 1
fi

sequences=()

for oeis in "$@"; do

	while [[ ${#oeis} -lt 6 ]]
		do oeis="0"$oeis
	done
	url="https://oeis.org/A$oeis/b$oeis.txt"

	# Download a sequence from OEIS
	echo "---"
	echo "Downloading the OEIS sequence A$oeis..."
	wget -N -q "$url"

	if [ $? -ne 0 ]; then
		echo "Can't download the file corresponding to A$oeis, skip it..."
		continue
	fi
	
	sequences+=("A$oeis")

	# Prepare the sequence as a file for factorization
	awk '{if ($1 ~ /^[0-9]+$/) print $2}' "b$oeis.txt" | sort -rn | uniq > "task.txt"

	# Process the "task.txt" file for factorization
	./factor --input-file "task.txt" --output-csv --verbose 2
	
	if [ $? -ne 0 ]; then
		# Exit with a message on error
		echo "Error during facotization of A$oeis"
		exit 1
	fi
done

echo "All OEIS sequences in [ "${sequences[*]}" ] were factored."
