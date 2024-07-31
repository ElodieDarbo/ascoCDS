#!/usr/bin/python

###########################################################
#                                                         #
#                      Compute GC%                        #
#                                                         #
###########################################################

# Syntax : python3 computeGC.py nuclGen.csv nuclCDS.csv

#----- Imports --------------------------------------------

import sys

#----- Main -----------------------------------------------

def main(argv):
	genCount = {}
	cdsCount = {}
	print( "\t".join(["#Species", "GCgen", "GCcod", "GCnc"]))
	
#----- Data input -----------------------------------------

	gen = open( argv[0], "r" )
	for line in gen:
		if not line.startswith("#"):
			line = line.rstrip()
			data = line.split("\t")
			genCount[data[0]] = [float(i) for i in data[1:]] # type conversion

	cds = open( argv[1], "r" )
	for line in cds:
		if not line.startswith("#"):
			line = line.rstrip()
			data = line.split("\t")
			cdsCount[data[0]] = [float(i) for i in data[1:]] # type conversion

#----- GC computations ------------------------------------

	for species in sorted(genCount.keys()):
		# GC% of the genome
		(Ag, Cg, Gg, Tg) = genCount[species]
		GCgen = (Cg + Gg) / (Ag + Cg + Gg + Tg)
		
		# GC% of the coding sequences
		(Ac, Cc, Gc, Tc) = cdsCount[species]
		GCcds = (Cc + Gc) / (Ac + Cc + Gc + Tc)
		
		# GC% of the non coding sequences
		GCnc = ((Cg + Gg) - (Cc + Gc)) / ((Ag + Cg + Gg + Tg) - (Ac + Cc + Gc + Tc))
		if GCnc > 1 or GCnc < 0:
			print(species, "Invalid GC", GCnc)
			exit()

#----- Output ---------------------------------------------

		print( "\t".join([species, str(GCgen), str(GCcds), str(GCnc)]))
		
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])