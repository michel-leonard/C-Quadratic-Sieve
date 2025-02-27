#!/bin/bash

# Ensure the factorization executable is available
if [[ ! -x "./factor" ]]; then
	echo "Error : The 'factor' executable isn't available."
	echo "Please place the executable in the current directory."
	exit 1
fi

# Define the OEIS sequences suitable for the factorization software
sequences=(
	'007262' #  1001 rows factored in  25 s
	'005575' #   898 rows factored in  25 s
	'007722' #    99 rows factored in  25 s
	'007257' #  1000 rows factored in  25 s
	'004982' #   401 rows factored in  24 s
	'006055' #    36 rows factored in  24 s
	'002559' #  1000 rows factored in  23 s
	'008280' #   872 rows factored in  23 s
	'004072' #     8 rows factored in  22 s
	'000889' #     7 rows factored in  22 s
	'005577' #   998 rows factored in  21 s
	'005709' #   495 rows factored in  21 s
	'000829' #     7 rows factored in  21 s
	'008421' #  1001 rows factored in  21 s
	'001531' #     8 rows factored in  20 s
	'002600' #  1000 rows factored in  20 s
	'001361' #  1000 rows factored in  20 s
	'008506' #  1001 rows factored in  20 s
	'008349' #  1001 rows factored in  19 s
	'001384' #   199 rows factored in  19 s
	'006289' #    22 rows factored in  19 s
	'006477' #   996 rows factored in  18 s
	'006271' #     6 rows factored in  18 s
	'007632' #   174 rows factored in  18 s
	'007252' #  1000 rows factored in  18 s
	'000856' #     7 rows factored in  17 s
	'002592' #   340 rows factored in  17 s
	'000877' #     7 rows factored in  17 s
	'000702' #   999 rows factored in  17 s
	'005973' #    29 rows factored in  16 s
	'003021' #   325 rows factored in  16 s
	'004071' #     8 rows factored in  16 s
	'000886' #     7 rows factored in  15 s
	'272919' #  1000 rows factored in  15 s
	'008392' #  1001 rows factored in  15 s
	'000610' #     9 rows factored in  15 s
	'001270' #  1248 rows factored in  15 s
	'000100' #   198 rows factored in  15 s
	'008362' #  1001 rows factored in  15 s
	'003260' #   197 rows factored in  15 s
	'003020' #   297 rows factored in  14 s
	'002590' #   274 rows factored in  14 s
	'008863' #  1001 rows factored in  14 s
	'001271' #  1141 rows factored in  14 s
	'008353' #   173 rows factored in  13 s
	'000545' #  1000 rows factored in  13 s
	'006733' #   158 rows factored in  13 s
	'002591' #   342 rows factored in  13 s
	'001537' #     5 rows factored in  13 s
	'000458' #   100 rows factored in  13 s
	'008377' #  1001 rows factored in  13 s
	'005710' #   494 rows factored in  13 s
	'007253' #   800 rows factored in  13 s
	'008394' #  1001 rows factored in  13 s
	'005422' #   298 rows factored in  12 s
	'008484' #   991 rows factored in  12 s
	'005036' #    99 rows factored in  12 s
	'000835' #     7 rows factored in  12 s
	'001590' #   198 rows factored in  12 s
	'006290' #    22 rows factored in  12 s
	'001467' #     7 rows factored in  12 s
	'008393' #  1001 rows factored in  12 s
	'008378' #  1001 rows factored in  11 s
	'002779' #   941 rows factored in  11 s
	'006053' #   198 rows factored in  11 s
	'000880' #     7 rows factored in  11 s
	'008376' #  1001 rows factored in  11 s
	'001005' #   198 rows factored in  10 s
	'008937' #   200 rows factored in  10 s
	'007263' #   667 rows factored in  10 s
	'008395' #  1001 rows factored in  10 s
	'006054' #   149 rows factored in  10 s
	'007654' #    99 rows factored in  10 s
	'008494' #  1001 rows factored in  10 s
	'001611' #   250 rows factored in   9 s
	'006732' #   158 rows factored in   9 s
	'006734' #   158 rows factored in   9 s
	'008504' #  1001 rows factored in   9 s
	'000832' #     7 rows factored in   9 s
	'006729' #   171 rows factored in   9 s
	'005199' #   118 rows factored in   9 s
	'007045' #   498 rows factored in   9 s
	'002574' #   197 rows factored in   9 s
	'002778' #   941 rows factored in   9 s
	'006424' #  1000 rows factored in   9 s
	'003819' #     6 rows factored in   8 s
	'000369' #  1225 rows factored in   8 s
	'007044' #   498 rows factored in   8 s
	'008295' #   493 rows factored in   8 s
	'007626' #  6863 rows factored in   8 s
	'000892' #     7 rows factored in   8 s
	'004534' #  5000 rows factored in   8 s
	'002721' #  1001 rows factored in   8 s
	'001328' #     8 rows factored in   8 s
	'006731' #   159 rows factored in   7 s
	'008505' #  1001 rows factored in   7 s
	'008922' #  1001 rows factored in   7 s
	'003458' #   374 rows factored in   7 s
	'008419' #  1001 rows factored in   7 s
	'008406' #   567 rows factored in   7 s
	'008417' #  1001 rows factored in   7 s
	'002599' #  1000 rows factored in   7 s
	'003821' #     5 rows factored in   7 s
	'000090' #   100 rows factored in   7 s
	'008398' #  1001 rows factored in   6 s
	'000518' #  1000 rows factored in   6 s
	'007258' #  1000 rows factored in   6 s
	'006730' #   171 rows factored in   6 s
	'006728' #   171 rows factored in   6 s
	'008495' #  1001 rows factored in   6 s
	'000585' #     8 rows factored in   6 s
	'000613' #     8 rows factored in   6 s
	'001326' #     8 rows factored in   6 s
	'006448' #     8 rows factored in   5 s
	'007994' #  1000 rows factored in   5 s
	'007827' #    98 rows factored in   5 s
	'000128' #   201 rows factored in   5 s
	'000126' #   201 rows factored in   5 s
	'004051' #  1000 rows factored in   4 s
	'006418' #   999 rows factored in   4 s
	'006727' #   172 rows factored in   4 s
	'001383' #   199 rows factored in   4 s
	'006429' #   999 rows factored in   4 s
	'000614' #     8 rows factored in   4 s
	'005922' #   200 rows factored in   4 s
	'001721' #   101 rows factored in   4 s
	'000235' #   197 rows factored in   4 s
	'006189' #   100 rows factored in   4 s
	'000814' #     7 rows factored in   4 s
	'000859' #     7 rows factored in   4 s
	'000618' #     9 rows factored in   4 s
	'002583' #   139 rows factored in   4 s
	'007571' #   148 rows factored in   3 s
	'007080' #    16 rows factored in   3 s
	'006675' #    99 rows factored in   3 s
	'003424' #  4738 rows factored in   3 s
	'002582' #   134 rows factored in   3 s
	'001578' #  1449 rows factored in   3 s
	'006691' #    28 rows factored in   3 s
	'000640' #   100 rows factored in   3 s
	'007505' #    27 rows factored in   3 s
	'004635' #  5000 rows factored in   3 s
	'000820' #     7 rows factored in   3 s
	'001716' #   101 rows factored in   3 s
	'003098' #   147 rows factored in   3 s
	'008608' #    26 rows factored in   3 s
	'002099' #   495 rows factored in   3 s
	'000339' #  7031 rows factored in   3 s
	'003824' #   516 rows factored in   3 s
	'004080' #   100 rows factored in   3 s
	'003508' #  2000 rows factored in   3 s
	'005516' #   988 rows factored in   3 s
	'002931' #    64 rows factored in   3 s
	'000409' #    13 rows factored in   3 s
	'002585' #    98 rows factored in   2 s
	'002584' #    99 rows factored in   2 s
	'000775' #   101 rows factored in   2 s
	'006737' #    81 rows factored in   2 s
	'006311' #    25 rows factored in   2 s
	'004174' #   698 rows factored in   2 s
	'007098' #    99 rows factored in   2 s
	'003381' #  8000 rows factored in   2 s
	'007546' #  5000 rows factored in   2 s
	'002185' #   569 rows factored in   2 s
	'002528' #    89 rows factored in   2 s
	'004636' #  5000 rows factored in   2 s
	'003149' #   101 rows factored in   2 s
	'007082' #    10 rows factored in   2 s
)

# Shuffle the OEIS sequences
sequences=( $(shuf -e "${sequences[@]}") )

for oeis in "${sequences[@]}"; do
	url="https://oeis.org/A$oeis/b$oeis.txt"

	# Download a sequence from OEIS
	echo "---"
	echo "Downloading the OEIS sequence A$oeis..."
	wget -N -q "$url"

	if [ $? -ne 0 ]; then
		echo "Can't download the file corresponding to A$oeis, skip it..."
		continue
	fi

	# Prepare the sequence as a file for factorization
	awk '{if ($1 ~ /^[0-9]+$/) print $2}' "b$oeis.txt" | sort -n | uniq > "task.txt"

	# Process the "task.txt" file for factorization
	./factor --input-file "task.txt" --output-csv --verbose 2
	
	if [[ $? -eq 0 ]]; then
		# Prompt the user to continue on success
		read -p "The sequence A$oeis is completed, press ENTER to factorize the next..."
		continue
	else
		# Exit with a message on error
		echo "Error during facotization of A$oeis"
		exit 1
	fi
done

echo "All OEIS sequences were processed."
