#!/usr/bin/env python

# Split a fasta or fastq file on a set of primer sequences.
# This script can be used to seperate multiplexed sequence data
# into seperate files.

# The script requires both an input sequence file (fasta or fastq)
# and a comma seperated primer file, containing the primer
# name and primer sequence.

# Arguments: Split_on_Primer.py -s [sequence file] -p [primer file] -r -m [mismatched nucleotides allowed] -c [number of threads]

# author: Youri Lammers
# contact: youri.lammers@naturalis.nl / youri.lammers@gmail.com

# import the modules used by the script
import os, argparse, itertools, sys

# Retrieve the commandline arguments
parser = argparse.ArgumentParser(description = 'Split a sequence file based on a list of primers.')

parser.add_argument('-f', '--sequence_file', metavar='Sequence file', dest='sequence', type=str,
			help='The sequence file in either fastq or fasta format.')
parser.add_argument('-p', '--primers', metavar='Primer file', dest='primer', type=str,
			help='Comma seperated value file containing the primers. Format = primer_name,primer_sequence')
parser.add_argument('-m', '--mis', metavar='Mismatches allowed', dest='mis', type=int,
			help='The maximum number of mismatches allowed between the primer and reads (default = 0)', default=0)
parser.add_argument('-s', '--shift', metavar='Nucleotide shift allowed', dest='shift', type=int,
			help='The maximum sequence shift allowed between the primer and reads (default = 0)', default=0)
parser.add_argument('-t', '--trim', action='store_true', dest='trim',
			help='Trim the primers from the sequences after splitting.')

args = parser.parse_args()

def extract_sequences():

	# open the sequence file submitted by the user
	sequence_file = open(args.sequence)

	# create a iterative index of all the headers
        lines = (x[1] for x in itertools.groupby(sequence_file, key=lambda line: line[0] == '@' or line[0] == '>'))

	# walk through the header and obtain the sequences (and quality score if applicable)
        for headers in lines:
                header = headers.next().strip()
                sequence = [line.strip() for line in lines.next()]

		# adjust for potential quality score lines starting with the @ symbol
                if len(sequence) == 2: sequence.append(headers.next().strip())

		# yield the header + sequence
		yield [header, sequence]


def extract_primers():

	# Split current primer list to a list and dictionary
	# List to keep track of # of hit per primer
	# Dictionary for easy access to output files

	# create the primer list, the list format is:
	# [primer_name, primer_squence+shift, primer_file, original_length]
	primer_list = []

	# set the output dictionary (same folders as the sequence file)
	# and get the extention for the ouput files (sames as input file)
	directory = os.path.dirname(os.path.realpath(args.sequence)) + '/'
	extension = os.path.splitext(args.sequence)[1]

	# walk through the primers in the primer file
	for line in open(args.primer):
		line = line.strip().split(',')

		# open the output_file
		output_file = open(directory+line[0]+extension,'w')
		
		# if the --shift argument > 0, create different primer
		# variants with the sequence shifts needed
		length = len(line[1])
		for i in range(0,args.shift+1):
			primer_list.append([line[0], line[1][i:], output_file, length])

	# create the unsorted file
	primer_list.append(['unsorted',0,open(directory+'unsorted'+extension,'w'),0])

	# return the primer dictionary
	return primer_list


def hamming_distance(sequence, primer):

	# calculate the hamming distance between the primer and read sequence
	# return the calculated distance
	return sum([s_nuc != p_nuc for s_nuc, p_nuc in zip(sequence, primer)])


def levenshtein_distance(sequence, primer):
	
	# calculate and return the levenshtein distance for the two sequences.
	# the levenshtein distance is only calculated if the hamming distance
	# was larger than the maximum mismatch and the --mis argument is
	# greater than 0 (default)

        previous = xrange(len(sequence) + 1)
        for pos_seq, nuc_seq in enumerate(sequence):
                current = [pos_seq + 1]
                for pos_prim, prim_seq in enumerate(primer):
                        insert, delete, change = previous[pos_prim + 1]+1, current[pos_prim]+1, previous[pos_prim] + (nuc_seq != prim_seq) 
                        current.append(min(insert, delete, change))
                previous = current

	# return the distance
        return previous[-1]


def trim_primer(sequence, length):

	# trim the primer from the sequence (either in fasta or fastq format)
	# function is only used if the --trim arguments is provided
	sequence[1][0] = sequence[1][0][length:]
	if len(sequence[1]) > 1:
		sequence[1][2] = sequence[1][2][length:]
	
	# return the trimmed sequence
	return sequence


def compare_sequences(sequence, primer_list, read_shift, method, distance_results, max_mis):

	# compare the sequence to the primers

	# parse through the primers in the primer dictionary
	for primer in primer_list:

		# skip the unsorted category
		if primer[0] == 'unsorted': continue
		
		# get the primer information
		primer_name, primer_sequence, primer_file, primer_length = primer
		primer_shift = primer_length - len(primer_sequence)

		# skip same shift comparisons between the read and primers unless the shift
		# equals zero, futhermore skip shift that have been carried out before for the
		# primer - read combination. (ie read shift 2, primer shift 1 equals 
		# read shift 1, primer shift 0)
		if read_shift >= 1 and (read_shift - primer_shift) < read_shift: continue

		# calculate the distance with either the hamming or levenshtein method		
		if method == 'hamming': 
			distance = hamming_distance(sequence[:len(primer_sequence)], primer_sequence) + abs(read_shift - primer_shift)
		else:
			distance = levenshtein_distance(sequence[:len(primer_sequence)], primer_sequence) + abs(read_shift - primer_shift)
		
		# append the distance results if they are lower than the mis threshold
		if distance <= max_mis: distance_results.append([distance, primer_file, primer_length])
		
		# stop calculating distances when the distance equals the shift
		if distance == 0: break
	
	# return the calculated distance_results list (sorted)
	return distance_results


def find_best_primer(read, primer_list):

	# worker thread for the distance calculations
	# this function will calculate the distance
	# between the pimers and sequence with either the
	# hamming or levenhstein method

	# list with the distance results
	distance_results = []

	# create the sequence for each potential read shift indicated by
	# the --shift arugment
	for read_shift in range(0,args.shift+1):

		# set the shifted sequence
		sequence = read[1][0][read_shift:]

		# obtain the distance results for all primers using
		# the hamming distance method
		distance_results = compare_sequences(sequence, primer_list, read_shift, 'hamming', distance_results, args.mis)

		# if the best result equals the current shift, stop further comparisons
		if len(distance_results) > 0:
			distance_results.sort()
			if distance_results[0][0] == 0: break

	# If no distance could be obtained with the hamming distance method, the Levenshtein
	# method is used. Method is only used if the maximum number of mismatches is larger than 1.
	if len(distance_results) == 0 and args.mis > 0:
		distance_results = compare_sequences(read[1][0], primer_list, 0, 'levenshtein', distance_results, args.mis)

	
	# pick the best match if there are multiple valid
	# primers found
	if len(distance_results) >= 1:

		# sort the distance results list (ascending based on distance)
		distance_results.sort()

		# pick the best primer
		primer = distance_results[0]

		# trim the read if --trim is selected
		if args.trim == True: read = trim_primer(read, primer[2])

		# save the read to the output file
		#write_read(read, primer[1])
	else:
		# no primer matches with the read, the read unsorted
		# will be written to the unsorted file
		#write_read(read, primer_list[-1][2])
		pass


def write_read(read, output_file):

	# write the read to the output_file in either fasta or fastq
	# format depending on the read
	if len(read[1]) > 1:
		output_file.write('{0}\n{1}\n+\n{2}\n'.format(read[0], read[1][0], read[1][2]))
	else:
		output_file.write('{0}\n{1}\n'.format(read[0], '\n'.join([read[1][0][i:i+60] for i in range(0, len(read[1][0]), 60)])))

	
def main ():

	# obtain the list with the primer names
	# sequences and mismatch shifts
	primer_list = extract_primers()

	# get the largest value in the primer list + max mismatch
	size = sorted(((value[3] + args.mis) for value in primer_list), reverse=True)[0]

	# parse through the sequences and find
	# the closest matching primer, write the results
	# to the output file for said primer.
	for read in extract_sequences():
		# if the read is shorter than the maximum primer size
		# write the read to the unsorted file
		if len(read[1][0]) <= size:
			write_read(read, primer_list[-1][2])
		else:
			find_best_primer(read, primer_list)
	

if __name__ == '__main__':
	main()

