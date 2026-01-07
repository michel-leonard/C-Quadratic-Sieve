# Usage: python3 script.py "demand.txt" "results.csv"

import csv
import math
import random

# Check for primality using the Miller-Rabin algorithm
def is_prime(n, k = 40):
	if n < 2 or (n != 2 and not n & 1):
		return False
	if n < 6:
		return True
	random_gen = random.SystemRandom()
	for _ in range(k):
		a = random_gen.randrange(2, n - 1)
		exp = n - 1
		while not exp & 1:
			exp >>= 1
		if pow(a, exp, n) == 1:
			continue
		while exp < n - 1:
			if pow(a, exp, n) == n - 1:
				break
			exp <<= 1
		else:
			return False
	return True


# Process the files
def process_files(input_files, output_files):
	input_numbers = [ ]
	
	# Read the input file
	with open(input_files, 'r') as f:
		for line in f.readlines():
			line = line.strip()
			if line and not line.startswith('#'):
				input_numbers.append(int(line.strip().split(' ')[0]))

	# Read the output file
	with open(output_files, 'r') as f:
		reader = csv.DictReader(f)
		results = {int(row['Input']): row['Factors'].split(';') for row in reader}

	total_lines = len(input_numbers)
	passed_count = 0
	failed_count = 0

	# Verify each number
	for num in input_numbers:
		if num in results:
			factors = results[num]
			all_factors_valid = True

			# Verify the factors
			for factor in factors:
				try:
					factor = int(factor)
					if factor != 1 and factor != -1 and not is_prime(abs(factor)):
						print(f"❌ The factor {factor} of {num} isn't prime.")
						all_factors_valid = False
				except ValueError:
					print(f"❌ Invalid value found for {factor} with the number {num}")
					all_factors_valid = False
			
			# Say "passed" or "failed"
			if all_factors_valid:
				# print(f"✅ {num} : Passed")
				passed_count += 1
			else:
				print(f"❌ {num} : Failed")
				failed_count += 1
		else:
			print(f"❌ The number {num} was not found in the output file {output_files}.")
			failed_count += 1

	# Summarize
	print(f"Verification completed: {total_lines} lines.")
	print(f"✅ {passed_count} success.")
	print(f"❌ {failed_count} errors.")
	if failed_count:
		exit(1)

# Main function (CLI)
if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Verify the prime factors")
	parser.add_argument("input_files", help="The input file containing the numbers to factor")
	parser.add_argument("output_files", help="The CSV output file containing the numbers factored")
	args = parser.parse_args()

	# Lancer le traitement
	process_files(args.input_files, args.output_files)
