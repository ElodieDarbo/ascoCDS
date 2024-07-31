#!/usr/bin/python

###########################################################
#                                                         #
#         Extracts preferred or avoided contexts          #
#                                                         #
###########################################################

# For each species, represented by a FASTA file of their 
# CDS, computes context related measures and report them.

# Moura G, Pinheiro M, Silva R, Miranda I, Afreixo V, 
# Dias G, Freitas A, Oliveira JL, Santos MA. Comparative
# context analysis of codon pairs on an ORFeome scale. 
# Genome Biol. 2005;6(3):R28. Epub 2005 Feb 15.

#----- Imports --------------------------------------------

import sys
import argparse
import os
import numpy as np
from Bio import SeqIO
from codons import CodonMap

#----- Functions ------------------------------------------
def loadMap(name, code):
	# load a table for codon conversion
	codon_map = CodonMap.CodonMap()
	return codon_map.get(name, code)

def residualAnalysis(labels, pairs, limit):
	# cf citation supra
	
	# margins
	byrow = pairs.sum(axis=1) * np.ones((len(labels), len(labels)))
	byrow = byrow.T
	bycol = pairs.sum(axis=0) * np.ones((len(labels), len(labels)))
	total = pairs.sum()
	
	# Pearson's statistics
	matR = (pairs - ((byrow * bycol) / total)) / np.sqrt((byrow * bycol) / total)
	
	# variance
	matV = (1 - (byrow / total)) * (1 - (bycol / total))
	
	# adjusted residuals
	matD = matR / np.sqrt(matV)
	
	# extraction of preferred or avoided contexts
	matP = np.nonzero(matD > limit)
	matA = np.nonzero(matD < -limit)
	
	prefIdx = list(zip(matP[0], matP[1]))
	avoiIdx = list(zip(matA[0], matA[1]))
	
	prefered = []
	for first, second in prefIdx:
		prefered.append(labels[first] + labels[second])
	avoided = []
	for first, second in avoiIdx:
		avoided.append(labels[first] + labels[second])
		
	return prefered, avoided

def tally(seq, pairs, labels):
	# counts contexts = pairs of consecutive codons

	k = 0
	while k < len(seq) - 3: # seq does not have stop codon
		pairs[labels.index(seq[k:k + 3]), labels.index(seq[k + 3:k + 6])] += 1
		k += 3 # shifts one codon at a time

	return

def wobblelize(seq, table):
	# transforms codons in wobbled codons
	wob = ''
	k = 0 # index
	while k < len(seq):
		wob = wob + table[seq[k:k + 3]]
		k += 3 # codons not overlapping trinucleotides
	return wob

#----- Main -----------------------------------------------

def main(argv):

#----- Command line parsing -------------------------------

	parser = argparse.ArgumentParser(description='Computes contexts from CDS.')
	parser.add_argument('-f', '--fasta', required=True, help='Name of the input FASTA file for CDS.')
	parser.add_argument('-o' ,'--outdir', required=True, help='Name of the output directory.')
	parser.add_argument('-s', '--species', required=True, help='Name of the species.')
	parser.add_argument('-g', '--gencode', required=True, help='Number of the genetic code.')
	args = parser.parse_args()

#----- Initialisations ------------------------------------

	completeTable = loadMap("complete", int(args.gencode))
	wobbleTable = loadMap("wobble", int(args.gencode))
	
	codonsC = sorted(set(list(completeTable.values())))
	#codonsW = sorted(set(list(wobbleTable.values()) + ['CUA','CUY','CUH'])) # every possible codon name
	codonsW = sorted(set(list(wobbleTable.values())))
	for stopCodon in ['UAG','UGA','UAA']:
		codonsC.remove(stopCodon)
		codonsW.remove(stopCodon)
		
	contextC = np.zeros((len(codonsC), len(codonsC))) # tables of contexts
	contextW = np.zeros((len(codonsW), len(codonsW)))
	
	threshold = 9 # signification threshold for residual analysis
		
#----- Context census -------------------------------------

	for cds in SeqIO.parse( args.fasta, "fasta" ): # all CDS
		# preliminary computations
		sequence  = str(cds.seq).upper() # on nucleotide sequence
		sequence  = sequence[:-3] # removes stop codon
		sequenceW = wobblelize(sequence, wobbleTable) # converts sequence
		
		# census
		tally(sequence, contextC, codonsC)
		tally(sequenceW, contextW, codonsW)
		
#----- Context exploitation -------------------------------

	preferC, avoidC = residualAnalysis(codonsC, contextC, threshold)
	preferW, avoidW = residualAnalysis(codonsW, contextW, threshold)
	
#----- Output ---------------------------------------------

	# build name of the output file
	(filepath, filename) = os.path.split(args.fasta)
	(shortname, extension) = os.path.splitext(filename)
	
	outname = shortname + "_countsC.csv"
	outFile = os.path.join(args.outdir, outname)
	head = "\t".join(codonsC)
	countTable = contextC.astype(int)
	np.savetxt(outFile, countTable, fmt='%6d', delimiter="\t", header=head, comments='#')
	
	outname = shortname + "_countsW.csv"
	outFile = os.path.join(args.outdir, outname)
	head = "\t".join(codonsW)
	countTable = contextW.astype(int)
	np.savetxt(outFile, countTable, fmt='%6d', delimiter="\t", header=head, comments='#')
	
	outname = shortname + "_contextC.csv"
	outFile = os.path.join(args.outdir, outname)
	out = open(outFile, "w")
	out.write("#Preferred contexts\n")
	output = "\n".join(preferC)
	out.write(output + "\n")
	out.write("#Avoided contexts\n")
	output = "\n".join(avoidC)
	out.write(output + "\n")
	out.close()
	
	outname = shortname + "_contextW.csv"
	outFile = os.path.join(args.outdir, outname)
	out = open(outFile, "w")
	out.write("#Preferred contexts\n")
	output = "\n".join(preferW)
	out.write(output + "\n")
	out.write("#Avoided contexts\n")
	output = "\n".join(avoidW)
	out.write(output + "\n")
	out.close()
	
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
