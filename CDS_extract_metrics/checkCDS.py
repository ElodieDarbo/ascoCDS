#!/usr/bin/python

###########################################################
#                                                         #
#                        Check CDS                        #
#                                                         #
###########################################################

#----- Initialisations ------------------------------------

import sys
import argparse
import os
import re

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

__version__ = "2.0"
indir = './'
outdir = './'

#----- Check routines -------------------------------------

def checkLength( entry ):
  if len( entry ) % 3 == 0:
    return ''
  else:
    return '3'

def checkBegin( entry ):
  if entry.seq.upper()[0:3] in ['ATG','CTG','TTG']:
    return ''
  else:
    return 'b'

def checkEnd( entry ):
  if entry.seq.upper()[-3:] in ['TAA','TAG','TGA']:
    return ''
  else:
    return 'e'

def checkInterrupt( entry ):
  for i in range(3, len( entry ) - 3, 3):
    if entry.seq.upper()[i:i + 3] in ['TAA','TAG','TGA']:
      return 'i'
  return ''

def checkAmbiguosity( entry ):
  for i in range(len( entry )):
    if entry.seq.upper()[i] not in ['A','C','G','T']:
      return 'a'
  return ''
           
#----- Main -----------------------------------------------

def main(argv):
  
#----- Command line arguments -----------------------------

  parser = argparse.ArgumentParser(description='Checks the validity of CDS in FASTA file.')
  parser.add_argument('-i', '--indir', help='Name of the input directory.')
  parser.add_argument('-o' ,'--outdir', help='Name of the output directory.')
  args = parser.parse_args()
  
#----- Output files ---------------------------------------

  inputPath = args.indir
  outputPath = args.outdir
  log = open('checkCDS_log.txt', 'w')
  cover = open('species_cover.csv', 'w')
  cover.write( "\t".join( [ "species", "raw_CDS", "filtered_CDS", "Cover\n" ] ))
  
#----- Scanning input directory ---------------------------
  
  for entry in os.listdir(inputPath):
    if entry.endswith('fasta'):
      species = os.path.splitext(entry)[0]
      sys.stderr.write( species + "\n" )
    
#----- Parsing FASTA files --------------------------------

      totalCDS = 0
      filteredCDS = 0
      
      outfile = open(args.outdir + entry, 'w') # FASTA output
      
      for record in SeqIO.parse( args.indir + entry, "fasta" ):
      
#----- CDS sequence checks --------------------------------

        checks = ''
        checks += checkLength( record ) # length is a multiple of 3
        checks += checkBegin( record ) # must begin with a start codon
        checks += checkEnd( record ) # must end with a stop codon
        checks += checkInterrupt( record ) # must not end with an in-phase stop codon
        checks += checkAmbiguosity( record) # must not contain ambiguous nucleotide
  
        if checks:
          output = "\t".join( [ record.id, species, checks ] ) + "\n"
          log.write( output ) # report rejection diagnosis

        else:
          RNA = SeqRecord( record.seq.transcribe() , id=record.id, description='' )
          # overriding SeqIO as it is not possible to write sequences one by one
          outfile.write(RNA.format("fasta")) 

          #SeqIO.write( RNA, args.outdir + entry, 'fasta' )
          filteredCDS += 1
          
        totalCDS += 1
      outfile.close()
      
#----- Species statistics ---------------------------------

      coverage = filteredCDS / totalCDS
      output = "\t".join( [ species, str(totalCDS), str(filteredCDS),  str(coverage) ] ) + "\n"
      cover.write(output)

#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
  

