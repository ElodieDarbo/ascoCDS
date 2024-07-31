#!/usr/bin/python

###########################################################
#                                                         #
#                    Count all codons                     #
#                                                         #
###########################################################

# For all species, represented by a FASTA file of their CDS,
# computes codon statistics

#----- Imports --------------------------------------------

import sys
import argparse
import os
from Bio import SeqIO
from Codon_map import codonmap

#----- Main -----------------------------------------------

def main(argv):

#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(description='Computes codon statistics.')
  parser.add_argument('-i', '--indir', help='Name of the input directory of FASTA files.')
  parser.add_argument('-o' ,'--outdir', help='Name of the output directory.')
  parser.add_argument('-s', '--species', help='Name of the table species / file (csv).')
  args = parser.parse_args()

#----- Initialisations ------------------------------------
  
  # builds a dictionary "species -> gencode"
  gencode = {}
  fh = open(args.species, "r")
  for line in fh:
    if not line.startswith("#"):
      line.rstrip()
      data = line.split("\t")
      gencode[data[0]] = int(data[1])
  fh.close()
  
  wobble = {}
  
#----- Species scan ---------------------------------------

  for species in sorted(gencode.keys()):
    comp_file = open(args.outdir + species + '.cod', 'w')
    wobb_file = open(args.outdir + species + '.wob', 'w')
  
    # make headers for the two output files, based on CodonMap
    #codon_map = CodonMap.CodonMap()
    codonCorresp = codonmap("wobble", gencode[species]) # dictionary with codons as keys and wobble as values
    
    codon_names = set(codonCorresp.keys())
    wobble_names = set(codonCorresp.values())
    for triplet in ['UAA','UAG','UGA']: # removes stop codons
      codon_names.discard(triplet)
      wobble_names.discard(triplet)
    for triplet in ['CUA', 'CUC', 'CUU', 'CUY', 'CUH']:
      wobble_names.add(triplet) # for gencodes 12 and 26
    
    # prints the headers in the result files
    comp_file.write("\t".join(["#CDS-id", "Species"] + list(sorted(codon_names))) + "\n")
    wobb_file.write("\t".join(["#CDS-id", "Species"] + list(sorted(wobble_names))) + "\n")

    print(species + "\t" + str(gencode[species]))
    sys.stdout.flush()
    
#----- Parsing sequences ----------------------------------
    
    for record in SeqIO.parse( args.indir + species + '.fasta', "fasta" ):
      sequence = str(record.seq).upper()

#----- Counting codons ------------------------------------

      complete_count = {}
      k = 0
      while k < len(sequence) - 3: # skips the stop codon
        if sequence[k:k+3] in complete_count.keys():
          complete_count[sequence[k:k+3]] += 1
        else:
          complete_count[sequence[k:k+3]] = 1
        k += 3 # next codon
        
      wobble_count = {}
      k = 0
      while k < len(sequence) - 3:
        if codonCorresp[sequence[k:k+3]] in wobble_count.keys():
          wobble_count[codonCorresp[sequence[k:k+3]]] += 1
        else:
          wobble_count[codonCorresp[sequence[k:k+3]]] = 1
        k += 3
        
#----- Computing codon frequencies ------------------------

      complete_freq = {}
      wobble_freq = {}
      total = ( len(sequence) - 3 ) / 3
      
      for codon in complete_count.keys():
        complete_freq[codon] = complete_count[codon] / total
        
      for codon in wobble_count.keys():
        wobble_freq[codon] = wobble_count[codon] / total    

#----- Table Completion -----------------------------------

      for codon in codonCorresp.keys(): # codons of the complete genetic code
        if codon not in complete_freq.keys():
          complete_freq[codon] = 0
      
      for codon in codonCorresp.values(): # codons of the wobbled genetic code
        if codon not in wobble_freq.keys():
          wobble_freq[codon] = 0
      for codon in ['CUA', 'CUC', 'CUU', 'CUY', 'CUH']:
        if codon not in wobble_freq.keys():
          wobble_freq[codon] = 0
          
#----- Output ---------------------------------------------

      data = [record.id, species]
      for codon in sorted(complete_freq.keys()):
        data.append(str(complete_freq[codon]))
      comp_file.write("\t".join(data) + "\n")
      
      data = [record.id, species]
      for codon in sorted(wobble_freq.keys()):
        data.append(str(wobble_freq[codon]))
      wobb_file.write("\t".join(data) + "\n")

#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
