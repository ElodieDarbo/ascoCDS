#!/usr/bin/python

# Syntax : python3 cure-doublons.py FastaFile > fixedFastaFile

import sys

count = {}
fh = open(sys.argv[1], "r")
for line in fh:
	if line.startswith(">"):
		line = line.rstrip()
		data = line.split(" ")
		if data[0] in count:
			count[ data[0] ] += 1
		else:
			count[ data[0] ] = 1

fh.close()

number = {}
for name in count.keys():
	number[name] = list( range( 1, count[name] + 1 ) )

fh = open(sys.argv[1], "r")
for line in fh:
	line = line.rstrip()
	if line.startswith(">"): # Defline
		data = line.split(" ")		
		if count[ data[0] ] == 1:
			name = data[0]
			desc = [name] + data[1:]
			print( " ".join( desc ) )
		else: # doublon
			num = number[data[0]].pop(0)
			name = data[0] + "-" + str(num)
			desc = [name] + data[1:]
			print( " ".join( desc ) )
	else:
		print(line)


fh.close()


