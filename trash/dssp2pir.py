#!/usr/bin/env python
from sys import argv
if len(argv)<2:
	print "Usage: ",argv[0], " file.dssp"
	exit()


h = open(argv[1],"r")

while 1:
	l=h.readline()
	if(l.startswith("  #  RESIDUE AA STRUCTURE")):
		break

seq = ""
ss = ""
for line in h.readlines():
	#data = line.split()
	if(line[13]!="!"):
		seq += line[13]
		ss += line[16].replace(" ","C").replace("B","C").replace("G","C").replace("I","C").replace("T","C").replace("S","C")
		

print "> ",argv[1],"\n\n",seq,"\n\n> SS\n\n",ss


