## Last modified: August 08, 2024 - SHA

# Read this to understand how to use several functions within Pyhton script:
# https://realpython.com/python-main-function/#use-if-__name__-__main__-to-control-the-execution-of-your-code

# Quick tutorial of pandas: https://pandas.pydata.org/docs/user_guide/10min.html

import os, sys
from Bio import SeqUtils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

import multiprocessing as mp

import gzip
from mimetypes import guess_type
from functools import partial

import numpy as np
import pandas as pd
from itertools import islice

# Usage: Calc_base_composition_stats_assembled_reads_v2.py <fastq file> "output_prefix_name"

# Reference AncSR1 and AncSR2 REBC1 nucleotide sequences, with unbiased degenerate codons for the RH library sites
sr1_fwd = Seq("TGTGCTGTTTGTTCTGATTATGCTTCTGGTTATCATTATGGNGTNTGGWSNTGYNNNNNNTGYNNNNNNTTCTTCAAGAGATCTATTCAAGGTCATAATGATTATATGTGTCCAGCT")
sr2_fwd = Seq("TGTTTAATTTGTGGTGATGAAGCTTCTGGTTGTCATTATGGTGTNCTNACNTGYNNNNNNTGYNNNNNNTTCTTCAAGAGAGCTGTTGAAGGTCAACATAATTATTTATGTGCTGGT")

sr1_rev = sr1_fwd.reverse_complement()
sr2_rev = sr2_fwd.reverse_complement()

# Complementary function to fix reads
def fix_assembled_read(record_info):
	# Rationale: We introduced a barcode and "frameshift" N's to increase diversity during sequencing --> Reads are stacked and positions are shifted.
	# Solution: Slice read sequence to match sr1 or sr2 sequence length (removes bin barcode and N's): nt_search identifies if the pattern (e.g. 'sr1_fwd') is within the query sequence ('read_seq'). 
	# It returns a list where the first item is the search pattern, followed by the position(s) of the matches.
	# If assembled read matched the reversed target sequence, return the reverse complement of the read (i.e., the forward sequence)

	record_id, record_seq, record_qual = record_info
	read_seq = Seq(record_seq)

	#Determine whether the assembled read contains the AncSR1 or AncSR2 sequence, and in which direction:
	direction_match = None
	if len(SeqUtils.nt_search(str(read_seq),str(sr1_fwd))) == 2 or len(SeqUtils.nt_search(str(read_seq),str(sr2_fwd))) == 2:
		direction_match = "forward"
	elif len(SeqUtils.nt_search(str(read_seq),str(sr1_rev))) == 2 or len(SeqUtils.nt_search(str(read_seq),str(sr2_rev))) == 2:
		direction_match = "reverse"

	r = Seq("")

	if direction_match == "forward":
		if len(SeqUtils.nt_search(str(read_seq),str(sr1_fwd))) == 2:
			start = SeqUtils.nt_search(str(read_seq),str(sr1_fwd))[1]
			end = start + len(sr1_fwd)
			r = read_seq[start:end]
		else:
			start = SeqUtils.nt_search(str(read_seq),str(sr2_fwd))[1]
			end = start + len(sr2_fwd)
			r = read_seq[start:end]

	elif direction_match == "reverse":
		if len(SeqUtils.nt_search(str(read_seq),str(sr1_rev))) == 2:
			start = SeqUtils.nt_search(str(read_seq),str(sr1_rev))[1]
			end = start + len(sr1_rev)
			r = read_seq[start:end].reverse_complement()
		else:
			start = SeqUtils.nt_search(str(read_seq),str(sr2_rev))[1]
			end = start + len(sr2_rev)
			r = read_seq[start:end].reverse_complement()

	return r

# Complementary function to compute the base composition stats
def base_comopsition_stats(fastq_file,outputFilename):
	cwd = os.getcwd() # get current working directory

	#Process the fastq file and determine whether it is gzip'd (processing gzip files requires less memory)
	encoding = guess_type(fastq_file)[1]  # determines the file extension
	if encoding is None:
		_open = open
	
	elif encoding == 'gzip':
		_open = partial(gzip.open, mode='rt')
	else:
		raise ValueError('Unknown file encoding: "{}"'.format(encoding))

	# Parallel processing of reads:
	#n_cores = mp.cpu_count() # max number of available CPUs
	n_cores = 10
	pool = mp.Pool(n_cores)

	array_reads = []
	counter = 0
	fastqs = FastqGeneralIterator(_open(fastq_file)) # Saves all records as a general iterator object (saves RAM and is faster)
	while counter < 11:
		n_reads = list(islice(fastqs, 100000)) # iterates by groups of 100,000 reads
		mapped_reads_n = pool.map(fix_assembled_read,n_reads,chunksize=1000) #Analyze batches of 1,000 reads in parallel
		mapped_reads_n = [i for i in mapped_reads_n if len(i) != 0] #remove unmapped reads
		array_reads += mapped_reads_n
		counter += 1

	# Convert mapped reads to data frame
	read_df = pd.DataFrame(array_reads)

	# Export Read array
	#filepath0 = os.path.join(cwd,outputFilename+'_readArray.csv')
	#read_df.to_csv(filepath0, index = False)
	print("Dimensions of read array: ",read_df.shape)
	#log_file.write(">Dimensions of read array: " + str(read_df.shape) + "\n") # Just as a sanity check

	# Build empty array to store the base composition per cycle (position of the read)
	base_comp = pd.DataFrame(np.zeros([read_df.shape[1],4],dtype='int'),columns=list("ACGT"))

	# Iterate over the columns (cycles) of the read and compute the frequency of each nucleotide
	for column in read_df:
		count_A = list(read_df.iloc[:,column]).count('A')
		count_C = list(read_df.iloc[:,column]).count('C')
		count_G = list(read_df.iloc[:,column]).count('G')
		count_T = list(read_df.iloc[:,column]).count('T')
		total = count_A + count_C + count_G + count_T

		freq_A = count_A / total
		freq_C = count_C / total
		freq_G = count_G / total
		freq_T = count_T / total

		counts = list(pd.Series([freq_A,freq_C,freq_G,freq_T]))

		# Replace cells in empty array with new frequencies (note that column in the 'read_array' is the row of the 'base_comp' array)
		base_comp.loc[column] = counts

	# save final base_comp array to a file
	filepath = os.path.join(cwd,outputFilename+'_base_composition.csv')
	base_comp.to_csv(filepath, index = False)	


# Main function
def main(fastq_file,outputFilename):
	base_comopsition_stats(fastq_file,outputFilename)

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])




