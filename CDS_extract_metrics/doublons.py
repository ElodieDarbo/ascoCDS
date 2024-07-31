#!/usr/bin/python
import sys

names = set()

fh = open(sys.argv[1], "r")

for line in fh:
	line = line.rstrip()
	if line in names:
	 	print(line)
	else:
		names.add(line)