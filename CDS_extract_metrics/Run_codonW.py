#!/usr/bin/python

###########################################################
#                                                         #
#                       Run codonW                        #
#                                                         #
###########################################################

#----- Imports --------------------------------------------

import sys
import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#----- Main -----------------------------------------------

def main(argv):

#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(description='Computes optimal codons on transcriptomes.')
  parser.add_argument('-i', '--infile', required=True, help='Name of the input FASTA file.')
  parser.add_argument('-o' ,'--outdir', required=True, help='Name of the output directory.')
  parser.add_argument('-s', '--species', required=True, help='Name of the species.')
  parser.add_argument('-g', '--gencode', required=True, help='Number of the genetic code.')
  parser.add_argument('-c', '--codonw', required=True, help='Path of codonW tool.')
  args = parser.parse_args()

#----- Initialisations ------------------------------------
  
  # Conversion of the genetic code to an internal option number for codonW
  genetcode = 0
  if args.gencode == 1:  
    genetcode = 0
  elif args.gencode == 12:
    genetcode = 8
  elif args.gencode == 26:  
    genetcode = 9
  
#----- Output files ---------------------------------------

  if not os.path.exists(os.path.join(args.outdir, 'summaries')):
    os.makedirs(os.path.join(args.outdir, 'summaries'))

  sandboxFile = os.path.join(args.outdir, "current_transcriptome")
 
#----- Sandbox setup --------------------------------------

  numcds = {}
  counter = 0
  records = []
    
  for record in SeqIO.parse(args.infile, "fasta" ): # for each CDS

    # converts record ID into a number to avoid truncation by CodonW
    counter += 1
    numcds[str(counter)] = record.id
          
    # rewrites U in T as asked by CodonW
    modseq = Seq(str(record.seq).upper()).back_transcribe()
    
    # puts all records together
    records.append(SeqRecord(modseq, id=str(counter), description=""))
  
  # outputs the modified sequences into a sandbox (current directory)
  SeqIO.write(records, sandboxFile + ".fasta", "fasta")
   
#----- CodonW computations --------------------------------
  
  # quality check
  codonw = subprocess.run([os.path.join(args.codonw, "codonw"), 
                          sandboxFile + ".fasta", 
                          "-code", str(genetcode),
                          "-nomenu", 
                          "-silent"], 
                          check=True)

  # correspondance analysis
  codonw = subprocess.run([os.path.join(args.codonw, "codonw"), 
                          sandboxFile + ".fasta", 
                          "-coa_cu",
                          "-code", str(genetcode),                            
                          "-nomenu", 
                          "-silent"], 
                          check=True) 
                            
  # optimal codons                        
  codonw = subprocess.run([os.path.join(args.codonw, "codonw"), 
                          sandboxFile + ".fasta", 
                          "-cai_file", "cai.coa",
                          "-cbi_file", "cbi.coa",
                          "-fop_file", "fop.coa",
                          "-all_indices",
                          "-code", str(genetcode),                            
                          "-nomenu", 
                          "-silent"], 
                          check=True) 

#----- results output -------------------------------------
    
  # saves the summary file computed by CodonW
  os.rename("summary.coa", os.path.join(args.outdir, 'summaries', args.species + "_summary.coa"))
    
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])

  # python3 Run_codonW.py -i test/In/ -o test/Out/ -s Species.csv -c /home/durrens/Software/codonW