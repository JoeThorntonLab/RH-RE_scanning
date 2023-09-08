# -*- coding: utf-8 -*-
import sys, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils
from Bio import Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import multiprocessing as mp
from functools import partial

import gzip
from mimetypes import guess_type

from itertools import islice
import numpy as np

#################################################
# This script is used for demultiplexing of NGS reads based on two components: 1) background library barcode (REBC), and 2) bin barcode sequences (BBC).
# The script will produce a file containing the reads that match *exactly* the REBC, and have *at most* one mismatch in the BBC. It will record and write the reads in the forward direction.
# The name of the file will be '<Background_library-barcode_bin-barcode>.fastq'
# The script sorts into X number of files as it process the fastq file, and processes the reads in parallel by chunks.
#
# Example of usage: python demultiplex_spikein_binBC_v9.py assembled_extraction_1000seqs.fastq bin_barcodes_forward_R1.txt spike_in_refSeq_forward.txt dmx_files.txt
#
# 'assembled_extraction_1000seqs.fastq' is the fastq file that is to be demultiplxed
# 'bin_barcodes_forward_R1.txt' is the file containing the bin-barcodes
# 'library_barcodes_forward.txt' is the file containing the library-barcodes
# 'dmx_files.txt' is a file containing the desired library+binBarcode combo files to be demultiplexed (each combo goes in one line)
#
# Format of barcode files: Each barcode goes in a separate line, and the barcode and sequence are separated by a tab
# The files may contain comments at the beginning using '#'
#
# Examples:
# spike-in:
# AncSR1_ERE	TTATGGTGTTTGGTCTTGTGAAGGTTGTAAAGCT
#
# Bin barcode:
# BB3	TGACTG
#################################################

# Amplified AncSR1 and AncSR2 nucleotide sequences of the Spike-in refSeqs, with degenerate nucleotides for RH positions (without primer anealing sites)
sr1_fwd = Seq("GCTTCTGGTTATCATTATGGNGTNTGGWSNTGYNNNNNBTGYNNVNNBTTCTTCAAGAGATCTATTCAAGGTCA")
sr2_fwd = Seq("GCTTCTGGTTGTCATYATGGTGTNYTNACNTGYNNNNNBTGYNNVNNBTTCTTCAAGAGAGCTGTTGAAGGTCA")

sr1_rev = sr1_fwd.reverse_complement()
sr2_rev = sr2_fwd.reverse_complement()

cwd = os.getcwd() # get current working directory

# Library amplification annealing sequences:
#SR1_P5_F/SR1_P7_F	TGTGCTGTTTGTTCTGATTAT (length = 21)
#
#SR1_P7_R/SR1_P5/R	TAATGATTATATGTGTCCAGCT (length = 22)
#
#
#
#SR2_P5_F/SR2_P7_F	TGTTTAATTTGTGGTGATGAA (length = 21)
#
#SR2_P7_R/SR2_P5_R	ACATAATTATTTATGTGCTGGT (length = 22)


# Complementary function to fix reads
def fix_assembled_read(read_seq,direction_match):
	# Rationale: We introduced a barcode and "frameshift" N's to increase diversity during sequencing --> Reads are stacked and positions are shifted.
	# Solution: Slice read sequence to match sr1 or sr2 sequence length (removes bin barcode and N's): nt_search identifies if the pattern (e.g. 'sr1_fwd') is within the query sequence ('read_seq'). 
	# It returns a list where the first item is the search pattern, followed by the position(s) of the matches.
	# If assembled read matched the reversed target sequence, return the reverse complement of the read (i.e., the forward sequence)
	# Capture reads based only on the amplified region (to deal with cross-amplification)
	# For each read, capture the 6 nucleotides introduced with the P5 primer + extend the sequence to cover the annealing sequences of the F/R library amp. primers (important to check for cross amplification of library primers)
	# Return the read in forward direction + the bin barcode + F/R seqs.

	r = Seq("")
	P5_primer_length = None
	P7_primer_length = None

	if direction_match == "forward":
		P5_primer_length = 21
		P7_primer_length = 22
		if len(SeqUtils.nt_search(str(read_seq),str(sr1_fwd))) == 2:
			start = SeqUtils.nt_search(str(read_seq),str(sr1_fwd))[1]
			end = start + len(sr1_fwd)
			r = read_seq[start-P5_primer_length-6 : end+P7_primer_length]
		else:
			start = SeqUtils.nt_search(str(read_seq),str(sr2_fwd))[1]
			end = start + len(sr2_fwd)
			r = read_seq[start-P5_primer_length-6 : end+P7_primer_length]

	elif direction_match == "reverse":
		P5_primer_length = 22
		P7_primer_length = 21
		# For reads matching in the rev direction, return reverse complement of read but keep bc sequence in the fwd direction
		if len(SeqUtils.nt_search(str(read_seq),str(sr1_rev))) == 2:
			start = SeqUtils.nt_search(str(read_seq),str(sr1_rev))[1]
			end = start + len(sr1_rev)
			bc = read_seq[start-P5_primer_length-6 : start-P5_primer_length]
			r = read_seq[start-P5_primer_length : end+P7_primer_length].reverse_complement()
			r = bc + r
		else:
			start = SeqUtils.nt_search(str(read_seq),str(sr2_rev))[1]
			end = start + len(sr2_rev)
			bc = read_seq[start-P5_primer_length-6 : start-P5_primer_length]
			r = read_seq[start-P5_primer_length : end+P7_primer_length].reverse_complement()
			r = bc + r

	return r

# Complementary function to return the Hamming distance between string1 and string2.
# string1 and string2 should be the same length.
def hamming_distance(string1, string2): 
	# Start with a distance of zero, and count up
	distance = 0
	# Loop over the indices of the string
	L = len(string1)
	for i in range(L):
		# Add 1 to the distance if these two characters are not equal
		if string1[i] != string2[i]:
			distance += 1
	# Return the final count of differences
	return distance


# Complementary function to assign reads to a specific bin barcode group
def assign_binBC(read_seq,barcode_file):
	#Extract bin barcode sequence from read:
	bc_in_read = read_seq[0:6]
	bc_assigned = None

	# Process barcode file
	with open(barcode_file,"r") as bc:
		for line in bc:
			line=line.strip()
			if not line.startswith("#"):
				data=line.split("\t")
				bin_barcode_name=data[0]
				bin_barcode_seq=data[1]
				
				# Allow at most 1 mismatch in the bin barcode
				aln_score = 0
				if len(str(bc_in_read)) > 0 and len(bin_barcode_seq) > 0:
					aln_score = 6 - hamming_distance(str(bc_in_read),bin_barcode_seq)
				# Perfect match should return a score of 6.
				if aln_score >= 5:
					bc_assigned = bin_barcode_name

	return bc_assigned


# Complementary function to assign reads to a specific background and RE library barcode
def assign_REBC(read_seq, REbarcode_file):

	background_assigned = None
	lib_bc_assigned = None

	# Process barcode file
	with open(REbarcode_file,"r") as bc:
		for line in bc:
			line=line.strip()
			if not line.startswith("#"):
				data=line.split("\t")
				lib_barcode_name=data[0]
				lib_barcode_seq=data[1]
				
				if len(SeqUtils.nt_search(str(read_seq),str(lib_barcode_seq))) == 2:
					lib_bc_assigned = lib_barcode_name

	result = str(lib_bc_assigned)
	return result


#Complemetary function to open sorting files: Create and open the files to sort the reads.
def open_files(dmx_files):
	dmx_file_list = []
	outfile_list = []
	with open(dmx_files, "r") as dmxf:
		for line in dmxf:
			line=line.strip()
			Rd_out = os.path.join(cwd,line + '.fasta')
			outfile_list.append(open(Rd_out, 'w')) # open new files
			dmx_file_list.append(line)

	return(dmx_file_list,outfile_list)


#Complementary function to map a read to a demultiplex file (on the fwd direction) based on REBC, BBC and 
# direction of match to reference (AncSR1 or AncSR2)
def map_read(record_info,REbarcode_file,barcode_file,dmx_file_list):
	#Extract the record info from the element stored in the generator as a tuple
	record_id, record_seq, record_qual = record_info
	read_name = record_id
	read_seq = Seq(record_seq)
	
	#Determine whether the assembled read contains the AncSR1 or AncSR2 sequence, and in which direction:
	direction_match = None
	if len(SeqUtils.nt_search(str(read_seq),str(sr1_fwd))) == 2 or len(SeqUtils.nt_search(str(read_seq),str(sr2_fwd))) == 2:
		direction_match = "forward"
	elif len(SeqUtils.nt_search(str(read_seq),str(sr1_rev))) == 2 or len(SeqUtils.nt_search(str(read_seq),str(sr2_rev))) == 2:
		direction_match = "reverse"
	
	#fix read and save as new sequence record --> get read in forward direction, with barcode, and without N's:
	REBC_info = None
	binBC_info = None
	mapped_record = None
	index_file = None
	if direction_match is not None:
		fixed_read = fix_assembled_read(read_seq,direction_match) 
	
		#get barcode information of the read:
		REBC_info = assign_REBC(fixed_read,REbarcode_file)
		binBC_info = assign_binBC(fixed_read,barcode_file)
		read_assignation = str(REBC_info) + '_' + str(binBC_info)

		#assign read to corresponding file
		index = [dmx_file_list.index(i) for i in dmx_file_list if read_assignation in i] #in which file should the read go
		# If read info matches a file name in 'dmx_files', save to that file
		if len(index) != 0 and len(index) == 1:
			mapped_record = SeqRecord(fixed_read, id=read_name, description = "")
			index_file = index[0]

	return [mapped_record, index_file]


# Complementary function to sort reads into a specific file
def sort_reads(file_index,outfile_list,mapped_reads):
	reads_to_sort = [i[0] for i in mapped_reads if i[1] == file_index] #filter reads that match index file
	SeqIO.write(reads_to_sort, outfile_list[file_index], "fasta")


# Main function to perform dmx on a file. The main argument is a file containing the 'REBC_BBC's combos ('dmx_files').
def main(fastq_file,barcode_file,REbarcode_file,dmx_files):
	
	#Process the fastq file and determine whether it is gzip'd (processing gzip files requires less memory)
	encoding = guess_type(fastq_file)[1]  # determines the file extension
	if encoding is None:
		_open = open
	
	elif encoding == 'gzip':
		_open = partial(gzip.open, mode='rt')
	else:
		raise ValueError('Unknown file encoding: "{}"'.format(encoding))

	#Create and open the files to sort the reads.	
	files = open_files(dmx_files)
	dmx_file_list = files[0]
	outfile_list = files[1]

	#set constant values to all arguments of 'map_read' function which are not changed during parallel processing
	partial_map_read = partial(map_read,REbarcode_file=REbarcode_file,barcode_file=barcode_file,dmx_file_list=dmx_file_list)

	# Parallel processing of reads:
	#n_cores = mp.cpu_count() # max number of available CPUs
	n_cores = 10
	pool = mp.Pool(n_cores)

	indexed_reads = []
	fastqs = FastqGeneralIterator(_open(fastq_file)) # Saves all records as a general iterator object (saves RAM and is faster)
	while True:
		n_reads = list(islice(fastqs, 100000)) # iterates by groups of 100,000 reads
		if not n_reads:
			break
		mapped_reads_n = pool.map(partial_map_read,n_reads,chunksize=1000) #Analyze batches of 1,000 reads in parallel
		mapped_reads_n = [i for i in mapped_reads_n if i[1] != None] #remove unmapped reads
		indexed_reads += [i[1] for i in mapped_reads_n] # indexes of files assigned to reads

		#sort reads into each file
		for i in [dmx_file_list.index(i) for i in dmx_file_list]:
			sort_reads(i,outfile_list,mapped_reads_n)

	for file in outfile_list:
		file.close()
	pool.close()

	#Count number of reads assigned to each file
	file_indexes = [dmx_file_list.index(i) for i in dmx_file_list]
	reads_per_file = [indexed_reads.count(i) for i in file_indexes]

	dmx_file_list = np.array(dmx_file_list)
	reads_per_file = np.array(reads_per_file)
	tbl=np.column_stack((dmx_file_list,reads_per_file))
	tbl=tbl.tolist()

	print(">Reads saved to: ")
	print('\n'.join(['\t'.join(['{:4}'.format(item) for item in row]) for row in tbl]))


if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
