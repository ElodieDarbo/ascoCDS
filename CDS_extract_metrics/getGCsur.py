#!/usr/bin/python

###########################################################
#                                                         #
#                       getGCsur.py                       #
#                                                         #
###########################################################

# Verification if all CDS are recognised
# Syntax : python3 getGCsur.py V coordinateFile filteredCDSfile prefix

# Extraction of surrounding sequences and GC computation
# Syntax : python3 getGCsur.py E coordinateFile filteredCDSfile genomeFile

#----- Imports --------------------------------------------

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

#----- Global values --------------------------------------

distance = 50
coor = {}
cds = set()

#----- Main -----------------------------------------------

def main(argv):
  
#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(description='Computes GC% or surrounding sequence for every CDS.')
  parser.add_argument('-c', '--coord', help='Name of the CDS coordinate file.')
  parser.add_argument('-f' ,'--fCDS', help='Name of the fiterded CDS file.')
  parser.add_argument('-g', '--genome', help='Name of the genome FASTA file.')
  parser.add_argument('-p', '--prefix', help='Prefix for CDS names.')
  parser.add_argument('-t', '--treatment', help='Treatment : [V]erify or [E]xtract.')
  args = parser.parse_args()
  
  if not args.coord:
    print("***** NAME OF COORDINATE FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  if not args.fCDS:
    print("***** NAME OF FILTERED CDS FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  if not args.treatment:
    print("***** CHOOSE TREATMENT *****")
    parser.print_help()
    sys.exit()
  if args.treatment == "E" and not args.genome:
    print("***** NAME OF GENOME FASTA FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  if not args.prefix:
    args.prefix = ''
  
#----- Coordinates ----------------------------------------

  coo = open(args.coord, "r")
  for line in coo:
    if not line.startswith("#"):
      line = line.rstrip()
      data = line.split("\t")
      coor[args.prefix + data[0]] = data[1:]
  coo.close()
  
#----- Verification ---------------------------------------

  if args.treatment == "V":
    for record in SeqIO.parse( args.fCDS, "fasta" ):    
      if record.id not in coor.keys():
        print(record.id, "not found in ", args.coord)
        exit()
        
#----- Filtered CDS ---------------------------------------

  if args.treatment == "E":
    for record in SeqIO.parse( args.fCDS, "fasta" ):
      if record.id not in coor.keys():
        print(record.id, "not found in ", args.coord)
        exit()
      else:
        cds.add(record.id) # set of all filtered CDS
        
#----- Extraction -----------------------------------------

    print("\t".join(["#CDSid", "GCsur"])) # header line
    for record in SeqIO.parse( args.genome, "fasta" ):
      
      # selection of the CDS present in the current chromosome
      relevantCDS = [ Id for Id in cds if coor[Id][0] == record.id ]
      
      for ID in relevantCDS:
        
        # boundaries for sequence extraction
        upborder = int( coor[ID][1] ) - distance - 1
        downborder = int( coor[ID][2] ) + distance

        # extraction of surrounding sequences
        if upborder >= 0:
          upstream = record[ upborder : int(coor[ID][1]) - 1 ].seq
        else:
          upstream = record[ : int(coor[ID][1]) - 1 ].seq

        if downborder < len( record ):
          downstream = record[ int(coor[ID][2]) : downborder ].seq
        else:
          downstream = record[ int(coor[ID][2]) : ].seq

				# GC% computation
        sur = upstream + downstream
        if len( sur ) >= distance:
          GCsur = str( ( sur.count("C") + sur.count("G") ) / len(sur) )
        else:
          GCsur = "NA"

				# output
        print("\t".join([ID, GCsur]))	

#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
