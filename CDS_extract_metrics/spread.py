#!/usr/bin/python

# Syntax : python3 spread.py -k key -o fromFile -t toFile

import sys

#----- Imports --------------------------------------------

import sys
import argparse

#----- Global values --------------------------------------


#----- Main -----------------------------------------------

def main(argv):
  
#----- Command line parsing -------------------------------

  parser = argparse.ArgumentParser(description='Spread species measures to sequence measures.')
  parser.add_argument('-k', '--key', help='Key to match data.')
  parser.add_argument('-o' ,'--origin', help='Name of the species measure file.')
  parser.add_argument('-t', '--to', help='Name of the sequence measure file.')
  args = parser.parse_args()
  
  if not args.origin:
    print("***** NAME OF SPECIES MEASURE FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  if not args.to:
    print("***** NAME OF SEQUENCE MEASURE FILE IS MISSING *****")
    parser.print_help()
    sys.exit()
  if not args.key:
    print("***** CHOOSE TREATMENT *****")
    parser.print_help()
    sys.exit()
  
#----- Species file ---------------------------------------

  head = ''
  addline = ''
  
  spf = open(args.origin, "r")
  for line in spf:
    line = line.rstrip()
    if line.startswith("#"):
      head = line.replace('#','')
    else:
      if args.key in line:
        addline = line
  spf.close()
  
#----- Sequence file --------------------------------------

  sef = open(args.to, "r")
  for line in sef:
    line = line.rstrip()
    if line.startswith("#"):
      print("\t".join([line, head]))
    else:
      print("\t".join([line, addline]))
      
#----- If called as script, call main program ------------- 

if __name__ == '__main__':
  main(sys.argv[1:])
