#!/usr/bin/python

###########################################################
#                                                         #
#                    getCOORDS_gff2.py                    #
#                                                         #
###########################################################

# Extraction of chromosome, start, end for all CDS from annotation file.
# Syntax : python3 getCOORDS_GFF2.py gff_file species acronym

#----- Imports --------------------------------------------

import sys

#----- Main -----------------------------------------------

def main(argv):
	GFFfile = argv[0]
	species = argv[1]
	acronym = argv[2]
	
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
						
#----- Output ---------------------------------------------

	for CDSid in CDSlines.keys():

		start = min([ x[3] for x in CDSlines[CDSid] ]) # start of the first exon
		end = max([ x[4] for x in CDSlines[CDSid] ]) # end of the last exon

		print("\t".join([CDSid, CDSlines[CDSid][0][0], start, end]))

#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])