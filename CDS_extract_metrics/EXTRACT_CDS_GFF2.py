#!/usr/bin/python

###########################################################
#                                                         #
#                    EXTRACT_CDS_gff2.py                  #
#                                                         #
###########################################################

# Extraction of all CDS from annotation file, output in FASTA format
# Syntax : python3 EXTRACT_CDS_GFF2.py gff_file fasta_file species acronym

#----- Imports --------------------------------------------

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

#----- Global values --------------------------------------


#----- Functions ------------------------------------------


#----- Main -----------------------------------------------

def main(argv):
	GFFfile = argv[0]
	FASTAfile = argv[1]
	species = argv[2]
	acronym = argv[3]
	
	sys.stderr.write( species + "\n" )
	
#----- Read CDS lines -------------------------------------

	CDSlines = {}
	gff = open( GFFfile, "r" )
	for line in gff:
		if not line.startswith( "#" ):
			line = line.rstrip()
			fields = line.split("\t")
			if fields[2] == "CDS": # line describing CDS
				attributes = fields[8].split("; ")
				for attribute in attributes:
					if attribute.startswith("proteinId"):
						CDSid = acronym + "_" + attribute.split(" ")[1]
						if CDSid not in CDSlines.keys():
							CDSlines[CDSid] = []
						CDSlines[CDSid].append(fields)
						
#----- CDS extraction -------------------------------------

	# chromosome by chromosome
	for record in SeqIO.parse( FASTAfile, "fasta" ):
		sys.stderr.write( "\t" + record.id + "\n" )
		relevantCDS = [ id for id in CDSlines.keys() if record.id in CDSlines[id][0] ]
		
		# then exon by exon in a CDS
		for CDSid in relevantCDS:
			exons = sorted(CDSlines[CDSid], key=lambda x:int(x[3])) # sort by start value
			
			CDS = Seq( '', generic_dna )
			for exon in exons:
				
				# concatenation of exons
				CDS += record[ ( int( exon[3] ) - 1): int( exon[4] ) ].upper()
				
			# reverse complement when needed
			if exons[0][6] == '-': # reverse strand
				CDS = CDS.reverse_complement()
				
#----- Output ---------------------------------------------

			CDS.id = CDSid
			CDS.name = CDSid
			CDS.description = ''
			
			print( CDS.format("fasta"), end = '' )							
		
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])