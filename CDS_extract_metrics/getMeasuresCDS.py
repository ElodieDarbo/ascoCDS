#!/usr/bin/python

###########################################################
#                                                         #
#                      CDS measures                       #
#                                                         #
###########################################################

# For a species, represented by a FASTA file of its 
# CDS, computes several measures and report them.
# Results will be appended to a previously available
# results file, CDSid by CDSid, and saved with the same
# filename in the output directory. Species name is a single
# string where space is replaced by _.

# CDS sequence dependent measures (GCcds, CpG, GC1, GC2, GC3, G3skew, A3skew)
# Syntax : python3 getMeasureCDS.py -f CDSfastaFile \
#                                   -r resultFile \
#                                   -o outputDirectory/ \
#                                   -s species
#                                   -t A

# Codon dependent measures (Lcod, F, Nc, Psyn, SCUO)
# Syntax : python3 getMeasureCDS.py -f CDSfastaFile \
#                                   -r resultFile \
#                                   -o outputDirectory/ \
#                                   -s species
#                                   -t B

# Codon dependant measures relying on CodonW analysis
# Syntax : python3 getMeasureCDS.py -f CDSfastaFile \
#                                   -r resultFile \
#                                   -o outputDirectory/ \
#                                   -s species
#                                   -c codonW_summary  
#                                   -t C

# Context dependant measures
# Syntax : python3 getMeasureCDS.py -f CDSfastaFile \
#                                   -r resultFile \
#                                   -o outputDirectory/ \
#                                   -s species
#                                   -k contextFile 
#                                   -t D

# Protein properties
# Syntax : python3 getMeasureCDS.py -f CDSfastaFile \
#                                   -r resultFile \
#                                   -o outputDirectory/ \
#                                   -s species
#                                   -t E


# ----- F [Stop risk factor] ------------------------------
# Schmid P, Flegel WA (2011) Codon usage in vertebrates is
# associated with a low risk of acquiring nonsense 
# mutations. Journal of translational medecine 9:87.

# ----- Nc [Effective number of codons] -------------------
# Xiaoyan S, Qun Y, Xuhua X (2012) An improved 
# implementation of effective number of codons (Nc).
# Mol. Biol. Evol. 30:191-196.
# -----
# Satapathy S S, Sahoo A K, Ray S K, Ghosh T C (2017)
# Codon degeneracy and amino acid abundance influence the
# measures of codon usage bias: improved Nc and ENCprime
# measures. Genes to Cells 22:277-283.
# -----
# Novembre J A (2002) Accounting for Background 
# Nucleotide Composition When Measuring Codon Usage Bias.
# Mol. Biol. Evol. 19(8):1390-1394.

# ----- CAI [Codon Adaptation Index] ----------------------
# Sharp P M, Li W H (1987) The codon adaptation index - a 
# measure of directional synonymous codon usage bias, and 
# its potential applications.
# Nucleic Acids Research 15: 1281-1295.

# ----- CBI [Codon Bias Index] ----------------------------
# Bennetzen J L, Hall B D (1982) Codon selection in yeast. 
# J. Biol. Chem. 257: 3026–3031.

# ----- SCUO [Synonymous Codon Usage Order] ---------------
# Wan X F, Xu D, Kleinhofs A, Zhou J (2004) Quantitative 
# relationship between synonymous codon usage bias and GC
# composition across unicellular genomes. BMC Evol. Biol.
# 4:19.

# ----- Gravy [Hydropathicity] ----------------------------
# Kyte J, Doolittle R F (1982) A simple method for 
# displaying the hydropathic character of a protein".
# Journal of Molecular Biology. 157: 105–132.

# ----- Aromo [frequency of aromatic amino acids] ---------
# Lobry J R, Gautier C (1994) Hydrophobicity, expressivity 
# and aromaticity are the major trends of amino-acid usage
# in 999 Escherichia coli chromosome-encoded genes.
# Nucleic Acids Research 22: 3174-3180

# ----- II [instability index] ----------------------------
# Guruprasad K, Reddy B V B, Pandit M W (1990) Correlation 
# between stability of a protein and its dipeptide 
# composition: a novel approach for predicting in vivo 
# stability of a protein from its primary sequence.
# Protein Engineering 4:155-161

#----- Imports --------------------------------------------

import sys
import argparse
from argparse import RawTextHelpFormatter
import os
import re
import math 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from codons import CodonMap

#----- Functions ------------------------------------------

def addresults(folder, results, head, measures, suffix, species, fasta):
  # build name of the output file
  (filepath, filename) = os.path.split(results)
  (shortname, extension) = os.path.splitext(filename)
  outname = shortname + suffix + extension
  outFile = os.path.join(folder, outname)
  
  # add the new results
  res = open(results, 'r') # already existing measures
  out = open(outFile, 'w')
     
  for line in res:
    line = line.rstrip()
    if line.startswith("#"): # header
      out.write( line + "\t" + head + "\n" ) # updated header
      
    else: # existing measures
      data = line.split("\t")
      if data[0] in measures.keys(): # first field in measure file must be CDSid
        if data[1] == species: # second field must be species name
          
          newMeasures = [] # adds new measures
          for measure in measures[data[0]]: # converts to strings if needed
            if type(measure) is str:
              newMeasures.append(measure)
            else:
              newMeasures.append(str(measure))
          data.extend(newMeasures)
          
          output = "\t".join(data) # updated measures
          out.write( output  + "\n" )
          del measures[data[0]] # for verification, see below
          
        else: # inconsitency in species name
          sys.stderr.write("Error on species name (" + species + ") in file " + results + "\n")
          sys.exit()
      else: # inconstitency in CDSid
        sys.stderr.write("Error on CDSid present in " + results + " and not in " + fasta + "\n")
        sys.exit()
  
  if len(measures) > 0: # untreated item still present
    sys.stderr.write("Error on CDSid present in " + fasta + " and not in " + results + "\n")
    sys.exit()
 
  out.close()
  res.close()

def cai(nb, fit, syn):
  # codon adaptation index
  # see Sharp P M, Li W H (1987), p.1286, equation 8
  
  fitSum = 0.0
  length = 0
  for fam in syn.keys():
    if len(syn[fam]) > 1: # several synonyms in the family
      for codon in syn[fam]:
        if codon in nb.keys():
          fitSum += nb[codon] * math.log(fit[codon])
          length += nb[codon]
          
  return math.exp( fitSum / length )

def codonCount(seq):
  # counts occurences of each codon
  count = {} # initialisation
  k = 0 # index

  while k < len(seq):
    codon = seq[k:k + 3]
    if codon not in count.keys():
      count[codon] = 0 # initilisation on the fly
    count[codon] += 1 # counts codon
    k += 3 # codons not overlapping trinucleotides
  return count

def contexts(seq, prefr, avoid, syn):
  # Measures of context usage
  
  # Census of the preferred/avoided contexts in the CDS
  pCtxts = {} # preferred contexts of the CDS
  aCtxts = {} # avoided contexts of the CDS

  k = 0
  while k < len(seq) - 3: # seq does not have stop codon
    context = seq[k:k+6]
    if context in prefr:
      if context not in pCtxts.keys(): # initialisation
        pCtxts[context] = 0
      pCtxts[context] += 1
    
    elif context in avoid:
      if context not in aCtxts.keys(): # initialisation
        aCtxts[context] = 0
      aCtxts[context] += 1
    
    k += 3 # shifts one codon at a time
  
  # Lists of synonyms
  synFam = {}
  for fam in syn:
    for codon in syn[fam]:
      synFam[codon] = syn[fam]
      
  # Expected number of contexts
  frPref = {} # frequency of expected preferred contexts for then CDS
  frAvoi = {} # frequency of avoided contexts for then CDS  
  
  for context in pCtxts.keys():
    syn5 = synFam[context[0:3]] # synonyms of the context first codon
    syn3 = synFam[context[3:6]] # synonyms of the context second codon
    if len(syn5) * len(syn3) > 1: # there are synonym contexts
      
      # cartesian product of the 2 lists of synonyms
      synContexts = [ cod5 + cod3 for cod5 in syn5 for cod3 in syn3 ]    
      
      inProd = 0 # nb of preferred contexts in the cartesian product
      for ctxt in synContexts:
        if ctxt in prefr:
          inProd += 1
      frPref[context] = inProd / len(synContexts) # frequency of preferred contexts in the family
    
  for context in aCtxts.keys():
    syn5 = synFam[context[0:3]] # synonyms of the context first codon
    syn3 = synFam[context[3:6]] # synonyms of the context second codon
    if len(syn5) * len(syn3) > 1: # there are synonym contexts
      
      # cartesian product of the 2 lists of synonyms
      synContexts = [ cod5 + cod3 for cod5 in syn5 for cod3 in syn3 ]
      
      inProd = 0 # nb of avoided contexts in the cartesian product
      for ctxt in synContexts:
        if ctxt in avoid:
          inProd += 1
      frAvoi[context] = inProd / len(synContexts) # frequency of avoided contexts in the family
    
  # Counts
  nbPr = 0 # number of preferred contexts
  nbAv = 0 # number of avoided contexts
  nbPrRan = 0 # number of expected preferred contexts
  nbAvRan = 0 # number of expected avoided contexts
  nbTot = (len(seq) / 3) - 1 # total number of contexts
  
  for context in pCtxts.keys():
    nbPr += pCtxts[context]
    if context in frPref.keys():
      nbPrRan += pCtxts[context] * frPref[context]
    
  for context in aCtxts.keys():  
    nbAv += aCtxts[context]
    if context in frAvoi.keys():
      nbAvRan += aCtxts[context] * frAvoi[context]

  # Measures
  fpc = nbPr / nbTot # frequency of preferred contexts
  fac = nbAv / nbTot # frequency of avoided contexts
  boc = fpc - fac # balance of frequencies
  
  ipc = ( nbPr - nbPrRan ) / ( nbTot - nbPrRan ) # index of preferred contexts
  iac = ( nbAv - nbAvRan ) / ( nbTot - nbAvRan ) # index of avoided contexts
  bic = ipc - iac # balance of indexes
  
  return fpc, fac, boc, ipc, iac, bic
  
def effective(seq, syn, count):
  # computes Nc : effective number of codons
  
  # first computes a table of normalized expected usage for each codon :
  # Normalized expected frequency of a codon is the product 
  # of the frequencies of the nucleotides in the codon 
  # divided by the sum of the expected frequencies of all 
  # the codons in a family.
  # Thus the sum of normalized expected frequencies in a family
  # equals 1.

  # nucleotide frequencies in the CDS
  freq = {}
  for nucl in "ACGU":
    freq[nucl] = seq.count(nucl) / len(seq)

  # codon expected frequencies for the CDS
  expected = {}
  for fam in syn.keys():
    total = 1e-10 # pseudocount to avoid division by zero
    for codon in syn[fam]:
      expected[codon] = 1.0 # initialisation
      for nucl in codon: # codon frequency
        if nucl == "Y": # in case of wobblelized sequence
          expected[codon] *= (freq["C"] + freq["U"])
        elif nucl == "H": # in case of wobblelized sequence
          expected[codon] *= (freq["A"] + freq["C"] + freq["U"])
        elif nucl == "N": # in case of wobblelized sequence
          expected[codon] *= (freq["A"] + freq["C"]  + freq["G"] + freq["U"])
        else:  
          expected[codon] *= freq[nucl]
      total += expected[codon]
    for codon in syn[fam]:
      expected[codon] /= total # normalization

  # computation of the effective number of codons
  Nc = 0.0
  for fam in syn.keys():
    total = 0.0
    for codon in syn[fam]: # counts for a family
      if codon in count.keys():
        total += count[codon]
      
    if total > 0.0: # the family is present in the CDS
      X2 = 0.0
      for codon in syn[fam]:
        if codon in count.keys():
          X2 += ((count[codon] / total) - expected[codon])**2 / expected[codon] # deviation
      usage = (X2 + 1) / len(syn[fam])
      Nc += 1 / usage
      
  return Nc

def eligible( seq):
  # checks whether seq contains all 4 nucleotides
  for nucl in "ACGU":
    if nucl not in seq:
      return False
  return True

def extractSummary(codfile):
  # extracts adapation values and optimal codons from 
  # summary file made by codonW on all the CDS of the
  # species
  
  # initialisations
  optValues = []
  optCodons = []
  adapt = {}
  
  # regular expressions for data extraction
  optMotif = r"^[123](,[123]){15}$"
  adaMotif = r"^([A-Z]){3} [A-Z]([A-Za-z]){2} "
  
  # data extraction
  summary = open(codfile, "r")
  for line in summary:
    line = line.rstrip("\t\n")
    if re.search(optMotif, line): # reads optimal codons
    # line example : 3,2,2,1,2,3,3,3,1,2,2,2,1,3,2,2\n
      linedata = line.split(",")
      optValues.extend(linedata)  
    
    elif re.search(adaMotif, line): # reads adaptation values
    # line example : GAG Glu 1257.0 1.0000000\tGGG Gly  223.0 0.1024345
      for part in line.split("\t"):
        partdata = part.split()
        adapt[partdata[0]] = float(partdata[2]) # Xi
  summary.close()
  
  # codons and index in the order set by CodonW  
  valNucl = {"U" : 1, "C" : 2, "A" : 3, "G" : 4}  
  for p1 in "UCAG": # codons in the order set by CodonW
    for p3 in "UCAG":
      for p2 in "UCAG":
        codon = p1 + p2 + p3
        codonIdx = (( valNucl[p1] - 1) * 16 ) + valNucl[p2] + (( valNucl[p3] - 1 ) * 4 )

        # extraction of optimal codons
        if optValues[ codonIdx - 1 ] == '3': # 3 means optimal codon, 1 is rare codon, 2 is other
          optCodons.append(codon)

  return adapt, optCodons      

def extractValue(name, datafile):
  idx = 0
  
  dat = open(datafile, 'r') # already existing measures
  line = dat.readline() # header line
  line = line.rstrip()
  data = line.split("\t")
  if name in data:
    idx = data.index(name) # find relevant index
    
    line = dat.readline() # first measure line
    line = line.rstrip()
    values = line.split("\t")
    found = values[idx] # get relevant value at index
    dat.close()
    return found
  
  else:
    print("Error name not found in header.")
    exit()

def fitness(adapt, syn):
  # computes fitness values of codons within synonym family
  # using codon counts from highly biases CDS
  
  # maximal occurence values within families
  maxCount = {}
  for fam in syn.keys():
    maxCount[fam] = 0.5 # avoids 0, value according to Sharp & Li
    for codon in syn[fam]:
      if adapt[codon] > maxCount[fam]:
        maxCount[fam] = adapt[codon]
  
  # codon fitness within families   
  fitness = {}
  for fam in syn.keys():
    for codon in syn[fam]:
      fitness[codon] = adapt[codon] / maxCount[fam]
  
  return fitness

def GC(seq):
  # computes GC% and CpG%
  GCcds = (seq.count("G") + seq.count("C")) / len(seq)      
  CpG = (seq.count("CG") * 2) / len(seq)
  return GCcds, CpG

def GCx(counts):
  # computes GC[123]% from codon counts
  
  # initialisations
  GC1 = 0.0
  GC2 = 0.0
  GC3 = 0.0
  total = 0
  
  for codon in counts.keys():
    if codon[0] in "CG":
      GC1 += counts[codon]
    if codon[1] in "CG":
      GC2 += counts[codon]
    if codon[2] in "CG":
      GC3 += counts[codon]
    total += counts[codon]
  return (GC1 / total, GC2 / total, GC3 / total)

def loadMap(name, code):
  # load a table for codon conversion
  codon_map = CodonMap.CodonMap()
  return codon_map.get(name, code)

def N3skew(counts):
  
  # computes G/C and A/U balances on coding strand
  
  # initialisations
  A3 = 0
  C3 = 0
  G3 = 0
  U3 = 0
  A3skew = 0.0
  G3skew = 0.0
  
  for codon in counts.keys():
    if codon[2] in "A":
      A3 += counts[codon]
    elif codon[2] in "C":
      C3 += counts[codon]
    elif codon[2] in "G":
      G3 += counts[codon]
    elif codon[2] in "TU":
      U3 += counts[codon]
    
    # avoiding division by zero  
    if G3 + C3 > 0:
      G3skew = (G3 - C3) / (G3 + C3)
    else:
      G3skew = 'NA'
      
    if A3 + U3 > 0:
      A3skew = (A3 - U3) / (A3 + U3)
    else:
      A3skew = 'NA'
      
  return G3skew, A3skew

def optiMeasures(nb, opt, syn):
  # Fop : frequency of optimal codons, CBI : codon bias index
  # computed only on codons belonging to families of synonyms with several members
  
  totCod = 0 # total count of CDS codons 
  totOpt = 0 # total count of CDS optimal codons 
  totRan = 0 # expected total count of CDS optimal codons
  
  for fam in syn.keys():
    if len(syn[fam]) > 1: # several synonyms in the family
      nCod = 0 # count of CDS codons belonging to the family
      nOpt = 0 # count of CDS optimal codons belonging to the familiy
      fOpt = 0 # amount of optimal codons in the family
      for codon in syn[fam]:
        if codon in nb.keys():
          nCod += nb[codon]
          if codon in opt:
            nOpt += nb[codon]
            fOpt += 1
      totCod += nCod
      totOpt += nOpt
      totRan += nCod * (fOpt / len(syn[fam])) # expected count of optimal codons
  
  Fop = totOpt / totCod
  CBI = ( totOpt - totRan ) / ( totCod - totRan )
  return Fop, CBI

def readContext(contextFile):
  # reads prefered and avoided contexts from a file
  prefr = set()
  avoid = set()
  bag = ''
  context = open(contextFile, "r")
  for line in context:
    line = line.strip()
    if line.startswith("#Preferred"):
      bag = 'prefr'
    if line.startswith("#Avoided"):
      bag = 'avoid'
    else:
      if bag == 'prefr' :
        prefr.add(line)
      else :
        avoid.add(line)
      
  if prefr.isdisjoint(avoid):
    return prefr, avoid
  else: # error interception, a context cannot be preferred and avoided 
    print("Context inconsistency")
    print(prefr.intersection(avoid))
    exit()

def readSuffix(resultFile):
  # determine if it's a C, W or A file
  (filepath, filename) = os.path.split(resultFile)
  (shortname, extension) = os.path.splitext(filename)
  if shortname.endswith('C'):
    return 'C'
  elif shortname.endswith('W'):
    return 'W'
  elif shortname.endswith('A'):
    return 'A'
  else:
    sys.stderr.write("Error on result file name : neither C, W nor A\n")
    sys.exit()

def riskFactor(counts, risk, synonyms):
  # computes F : stop risk factor

  # initialisations
  riskscore = 0.0
  randomriskscore = 0.0
  meanCount = {}

  # average number of codons assuming an unbiased usage of codons
  for fam in synonyms.keys(): # family of synonyms
    total = 0.0
    for codon in synonyms[fam]:
      if codon in counts.keys():
        total += counts[codon]
    meanCount[fam] = total / len(synonyms[fam])
  
  # stop risk computation
  for fam in synonyms.keys():
    for codon in synonyms[fam]:
      if codon in counts.keys():
        riskscore += (counts[codon] * risk[codon])
        randomriskscore += (meanCount[fam] * risk[codon])

  # output
  if randomriskscore == 0.0:
    return "NA"
  else:
    return riskscore / randomriskscore

def SCUO(counts, syn):
  # computes SCUO : non-randomness of codon usage, 
  #                 0 <= SCUO <= 1, 
  #                 SCUO = 0 means total randomness
  # using information theory

  # initialisations
  SCUO = 0.0
  Ntot = 0
  NcodInFam = {}

  # number of codons in relevant families of synonyms   
  for fam in syn.keys():
    if len(syn[fam]) > 1: # eligible family of synonyms
      NcodInFam[fam] = 0
      for codon in syn[fam]:
        if codon in counts.keys():
          NcodInFam[fam] += counts[codon] # number of codons in CDS which belong to the family
          Ntot += counts[codon] # total number of codons in eligible families
  
  # family enthropy
  for fam in NcodInFam.keys(): # eligible families
    H = 0.0 # initialisation of the family entropy
    for codon in syn[fam]: # all codons of the family
      if codon in counts.keys():
        p = counts[codon] / NcodInFam[fam] # frequency of the codon within its family
        H += -(p * math.log2(p)) # entropy of the codon
      
  # SCUO
    Hmax = (- math.log2(1 / len(syn[fam]))) # maximal entropy, i.e. uniform 
                                            # distribution of codons
    O = (Hmax - H) / Hmax # normalized information
    F = NcodInFam[fam] / Ntot # frequency of the family
    SCUO += F * O # SCUO of the family added to the global SCUO
  
  return SCUO

def synRate(counts, syn):
  # computes the number of codons of a CDS, which have at least one synonym
  nbSyn = 0
  for fam in syn.keys():
    if len(syn[fam]) > 1: # eligible family of synonyms
      for codon in syn[fam]:
        if codon in counts.keys():
          nbSyn += counts[codon]
  return nbSyn

def wobblelize(seq, table):
  # transforms codons in wobbled codons
  wob = ''
  k = 0 # index
  while k < len(seq):
    wob = wob + table[seq[k:k + 3]]
    k += 3 # codons not overlapping trinucleotides
  return wob

def wobConversion(adapt, optim, table):
  # converts codonW data according to wobble
  wobAdapt = {}
  wobOptim = set()
  
  for codon in adapt:
    if table[codon] not in wobAdapt.keys():
      wobAdapt[table[codon]] = 0.0 # initialisation
    wobAdapt[table[codon]] += adapt[codon] # merging values
  
  for codon in optim:
    wobOptim.add(table[codon]) # a set removes duplicates
    
  return wobAdapt, list(wobOptim)

#----- Main -----------------------------------------------

def main(argv):

#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description='Computes measures on CDS.')
  parser.add_argument('-f', '--fasta', required=True, help='Name of the input FASTA file for CDS.')
  parser.add_argument('-o' ,'--outdir', required=True, help='Name of the output directory.')
  parser.add_argument('-s', '--species', required=True, help='Name of the species.')
  parser.add_argument('-r', '--results', required=True, help='Name of the file for measures done so far.')
  parser.add_argument('-t', '--task', required=True, type=str.upper, help='Task to be done (letter).\n\tA : sequence dependent measures (GCcds, CpG, GC1, GC2, GC3, G3skew, A3skew)')
  parser.add_argument('-c', '--codonw', help='name of the summary file from CodonW.')
  parser.add_argument('-k', '--context', help='name of the file for preferred and avoided contexts.')
  args = parser.parse_args()

  if args.task == 'C' and not args.codonw:
    print("***** NAME OF CODONW SUMMARY FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
    
  if args.task == 'D' and not args.context:
    print("***** NAME OF CONTEXT FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  
#----- Initialisations ------------------------------------

  measuresC = {}
  measuresW = {}
  measuresA = {}
  genCode = int(extractValue("Gencode", args.results))
  
  wobbleTable = loadMap("wobble", genCode)
  synonymsC = loadMap("synonymcomplete", genCode) # Synonyms apply to all CDS of a species
  synonymsW = loadMap("synonymwobble", genCode)
  synonymsA = loadMap("synonymaa", genCode)

  if args.task == 'C': # data from codonW 
    adaptCount, optimals = extractSummary(args.codonw) # from codonW statistics
    adaptCountW, optimalsW = wobConversion(adaptCount, optimals, wobbleTable) # same for wobble
 
    fitness_c = fitness(adaptCount, synonymsC)
    fitness_w = fitness(adaptCountW, synonymsW)
    fitness_a = fitness(adaptCount, synonymsA)

  if args.task == 'D': # data from codonW
    preferredC, avoidedC = readContext(args.context) # contexts for the species
    preferredW, avoidedW = readContext(args.context.replace("contextC", "contextW"))
    
#----- Tasks ----------------------------------------------
  
  for cds in SeqIO.parse( args.fasta, "fasta" ): # all CDS
    # preliminary computations
    sequence  = str(cds.seq).upper() # on nucleotide sequence
    sequence  = sequence[:-3] # removes stop codon
    sequenceW = wobblelize(sequence, wobbleTable) # converts sequence

    nb_codons  = codonCount(sequence) # counts of codons
    nb_codonsW = codonCount(sequenceW)
    
    if args.task == 'A':  # CDS sequence dependent measures
      
      GCcds, CpG = GC(sequence) # GC% and CpG% in cds
      GC1, GC2, GC3 = GCx(nb_codons) # GC% at different positions in codons
      G3skew, A3skew = N3skew(nb_codons) # unbalances within single strand
      
      measuresC[cds.id] = [ GCcds, CpG, GC1, GC2, GC3, G3skew, A3skew ]

    if args.task == 'B':  # Codon dependent measures
      
      Lcod = len(sequence) / 3 # Length in codons
      
      riskMap = loadMap("stoprisk", genCode)
      F = riskFactor(nb_codons, riskMap, synonymsC) # Stop risk factor
      
      if eligible(sequence): # must contain all 4 nucleotides
        Nc_c = effective(sequence, synonymsC, nb_codons) # effective number of codons
        Nc_w = effective(sequence, synonymsW, nb_codonsW)
        Nc_a = effective(sequence, synonymsA, nb_codons)
      else:
        Nc_c = 'NA'
        Nc_w = 'NA'
        Nc_a = 'NA'
        
      Psyn_c = synRate(nb_codons, synonymsC) / Lcod # proportion of codon with synonyms
      Psyn_w = synRate(nb_codonsW, synonymsW) / Lcod
      Psyn_a = synRate(nb_codons, synonymsA) / Lcod

      SCUO_c = SCUO(nb_codons, synonymsC)
      SCUO_w = SCUO(nb_codonsW, synonymsW)
      SCUO_a = SCUO(nb_codons, synonymsA)

      measuresC[cds.id] = [ Lcod, F, Nc_c, Psyn_c, SCUO_c ]
      measuresW[cds.id] = [ Lcod, F, Nc_w, Psyn_w, SCUO_w ]
      measuresA[cds.id] = [ Lcod, F, Nc_a, Psyn_a, SCUO_a ]
    
    if args.task == 'C':  # Codon dependant measures relying on CodonW analysis
      
      # CAI : Codon Adaptation Index
      CAI_c = cai(nb_codons, fitness_c, synonymsC) 
      CAI_w = cai(nb_codonsW, fitness_w, synonymsW)
      CAI_a = cai(nb_codons, fitness_a, synonymsA)
      
      # Fop : frequency of optimal codons, CBI : codon bias index 
      Fop_c, CBI_c = optiMeasures(nb_codons, optimals, synonymsC)
      Fop_w, CBI_w = optiMeasures(nb_codonsW, optimalsW, synonymsW)
      Fop_a, CBI_a = optiMeasures(nb_codons, optimals, synonymsA)
            
      measuresC[cds.id] = [ CAI_c, CBI_c, Fop_c ]
      measuresW[cds.id] = [ CAI_w, CBI_w, Fop_w ]
      measuresA[cds.id] = [ CAI_a, CBI_a, Fop_a ]
    
    if args.task == 'D':  # Context dependant measures
      
      # frequencies of preferred/avoided contexts
      # FPC : Frequency of Preferred Contexts
      # FAC : Frequency of Avoided Contexts
      # BOC : Balance Of Context frequencies
      
      # indexes of preferred/avoided contexts
      # IPC : Index of Preferred Contexts (similar to CBI)
      # IAC : Index of Avoided Contexts
      # BIC : Balance Of Indexes
      
      fpcC, facC, bocC, ipcC, iacC, bicC = contexts(sequence, preferredC, avoidedC, synonymsC) # measures for the CDS
      fpcW, facW, bocW, ipcW, iacW, bicW = contexts(sequenceW, preferredW, avoidedW, synonymsW) # measures for the CDS
      fpcA, facA, bocA, ipcA, iacA, bicA = contexts(sequence, preferredC, avoidedC, synonymsA) # measures for the CDS
      
      measuresC[cds.id] = [ fpcC, facC, bocC, ipcC, iacC, bicC ]
      measuresW[cds.id] = [ fpcW, facW, bocW, ipcW, iacW, bicW ]
      measuresA[cds.id] = [ fpcA, facA, bocA, ipcA, iacA, bicA ]
    
    if args.task == 'E':  # Protein properties
      translation = cds.seq.translate(table=genCode, to_stop=True) # converts CDS into protein sequence
      protein = str(translation)
      
      # measures from protein product
      protmes = ProteinAnalysis(protein) # BioPython routines
      Gravy = protmes.gravy() # average hydropathicity
      Aromo = protmes.aromaticity() # frequency of aromatic amino acids
      pI    = protmes.isoelectric_point() # isoelectric point
      II    = protmes.instability_index() # instability index
      
      measuresC[cds.id] = [ Gravy, Aromo, pI, II ]
      measuresW[cds.id] = [ Gravy, Aromo, pI, II ]
      measuresA[cds.id] = [ Gravy, Aromo, pI, II ]

#----- Output ---------------------------------------------

  if args.task == 'A':
    extrahead = "\t".join(['GCcds', 'CpG', 'GC1', 'GC2', 'GC3', 'G3skew', 'A3skew'])
    addresults(args.outdir, args.results, extrahead, measuresC, '', args.species, args.fasta)
    
  if args.task == 'B':
    extrahead = "\t".join(['Lcod', 'F', 'Nc', 'Psyn', 'SCUO'])
    addresults(args.outdir, args.results, extrahead, measuresC, '_C', args.species, args.fasta)
    addresults(args.outdir, args.results, extrahead, measuresW, '_W', args.species, args.fasta)
    addresults(args.outdir, args.results, extrahead, measuresA, '_A', args.species, args.fasta)

  if args.task == 'C':
    extrahead = "\t".join(['CAI', 'CBI', 'Fop'])    
    baseName = args.results[:len(args.results) - 6]
    addresults(args.outdir, baseName + '_C.csv', extrahead, measuresC, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_W.csv', extrahead, measuresW, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_A.csv', extrahead, measuresA, '', args.species, args.fasta)
    
  if args.task == 'D':  # Context dependant measures
    extrahead = "\t".join(['FPC', 'FAC', 'BOC', 'IPC', 'IAC', 'BIC'])
    baseName = args.results[:len(args.results) - 6]
    addresults(args.outdir, baseName + '_C.csv', extrahead, measuresC, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_W.csv', extrahead, measuresW, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_A.csv', extrahead, measuresA, '', args.species, args.fasta)
      
  if args.task == 'E':  # Protein properties
    extrahead = "\t".join(['Gravy', 'Aromo', 'pI', 'II'])
    baseName = args.results[:len(args.results) - 6]
    addresults(args.outdir, baseName + '_C.csv', extrahead, measuresC, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_W.csv', extrahead, measuresW, '', args.species, args.fasta)
    addresults(args.outdir, baseName + '_A.csv', extrahead, measuresA, '', args.species, args.fasta)
    
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
