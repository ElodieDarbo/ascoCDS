#!/usr/bin/python

###########################################################
#                                                         #
#                   Nucleotide census                     #
#                                                         #
###########################################################

# Syntax: python3 nuclCensus.py -i file -f format -s species
#----- Imports --------------------------------------------

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#----- Main -----------------------------------------------

def main(argv):

#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(description='Counts nucleotides.')
  parser.add_argument('-i', '--infile', help='Name of the input file.')
  parser.add_argument('-f' ,'--format', help='Format of the input file.')
  parser.add_argument('-s', '--species', help='Species name.')
  args = parser.parse_args()

#----- Initialisations ------------------------------------

  (nbA, nbC, nbG, nbT) = (0, 0, 0, 0)
  
#----- Parsing sequences ----------------------------------

  for record in SeqIO.parse( args.infile, args.format ):
      sequence = str(record.seq).upper()
      
      nbA += sequence.count("A")
      nbC += sequence.count("C")
      nbG += sequence.count("G")
      nbT += ( sequence.count("T") + sequence.count("U") )

#----- Output --------------------------------------------

  line = "\t".join( (args.species, str(nbA), str(nbC), str(nbG), str(nbT) ) )
  print( line )

#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
