#!/usr/bin/python

###########################################################
#                                                         #
#                    CORRECT_CDS.py                       #
#                                                         #
###########################################################

# Lengthen the CDS sequence
# Syntax : python3 CORRECT_CDS.py cds_file genome_file acronym

#----- Imports --------------------------------------------

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#----- Global values --------------------------------------

CDS = []
find = re.compile('loc:(\S+)\(([+-])\)([\d-]+)') 

#----- Main -----------------------------------------------

def main(argv):
	CDSfile = argv[0]     # FASTA format
	genomeFile = argv[1]  # FASTA format
	acronym = argv[2]

#----- Read CDS file --------------------------------------

	for record in SeqIO.parse( CDSfile, "fasta" ):
		CDS.append(record)
		
#----- Read genome file -----------------------------------

	for chromosome in SeqIO.parse( genomeFile, "fasta" ):

#----- Extract supplementary sequence ---------------------

		for record in CDS:
			match = re.search(find, record.description)
			chrom, strand, coor = match.groups()
			if chrom == chromosome.id:
				begin, end = coor.split('-')
				if strand == '+':
					patch = chromosome.seq[ int(end) : int(end) + 2 ].upper()
					corrected = record + patch
					corrected.id = acronym + record.id
					corrected.description = re.sub( coor,
												    begin + "-" + str( int( end ) + 2 ),
													record.description )
					print(corrected.format("fasta"), end="")
				else:
					patch = chromosome.seq[ int(begin) - 3: int(begin) - 1 ].upper()
					patch = patch.reverse_complement()
					corrected = record + patch
					corrected.id = acronym + record.id
					corrected.description = re.sub( coor,
												    str( int( begin ) - 3 ) + "-" + end,
													record.description )
					print(corrected.format("fasta"), end="")
						
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
	main(sys.argv[1:])
		