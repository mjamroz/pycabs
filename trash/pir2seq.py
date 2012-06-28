#!/usr/bin/env python
import sys


def findpos(seq, pat):

    matches = []
    current_match = seq.find(pat)
    while current_match != -1:
       matches.append(current_match)
       current_match =seq.find(pat, current_match+1)
    return matches


dict ={'A': 'ALA', 'R': 'ARG','N': 'ASN','D': 'ASP','C': 'CYS','E': 'GLU','Q': 'GLN','G': 'GLY','H': 'HIS','I': 'ILE','L': 'LEU','K': 'LYS','M': 'MET','F': 'PHE','P': 'PRO','S': 'SER','T':'THR','W': 'TRP','Y': 'TYR', 'V': 'VAL'}

dict_seq = {'H': 2, 'E': 4, '-': 1, 'C': 1, 'T': 1, 'S': 1, 'B': 1, 'G': 1, 'I': 1}



def print_seq_format(sequence, second_str):
	counter = 0
	for aminoacid in sequence:
		counter = counter + 1
		#print counter,aminoacid,second_str[counter-1], len(sequence), len(second_str)
		print '%5d   %s    %1d' % (counter,dict[aminoacid],dict_seq[second_str[counter-1]])
	return


buff = ""
sequence = ""
second_str = ""

if (len(sys.argv)==2):
       
	#print dict


	f = open(sys.argv[1],"r")

	f.readline()
	f.readline()

	while 1:				# czyta sekwencje
		line = f.readline()
                #print line
		if line[0]!=">":
			sequence  = sequence + line
		else:
			break

	sequence = sequence.replace("\n","").replace("*","").replace(" ","")

	while 1:				# czyta strukture
		line = f.readline()
		second_str  = second_str + line
		if (len(line)==0):
			break

	second_str  = second_str.replace("\n","").replace("*","").replace(" ","")

	#print sequence
	#print second_str
	print_seq_format(sequence, second_str)


	f.close()

###########################################################################################

elif (len(sys.argv)>2):
	structures = []
        sequences = []

	for i in range(1,len(sys.argv)):
		sequence = ""
		buff = ""

		f = open(sys.argv[i],"r")
		f.readline()
		f.readline()
		second_str = ""
		while 1:				# czyta sekwencje
			line = f.readline()
			#print line
			if line[0]!=">":
				sequence  = sequence + line
			else:
				break

			sequence = sequence.replace("\n","").replace("*","").replace(" ","")
                sequences.append(sequence)


		while 1:				# czyta strukture
			line = f.readline()
			second_str  = second_str + line
			if (len(line)==0):
				break
			second_str  = second_str.replace("\n","").replace("*","").replace(" ","")
		structures.append(second_str)
		f.close()
	
	control=sequences[0]
	for s in sequences[1:]:
	   #print s
	   if s!=control:
	      print "Oops! sequence lengths are not equal, bye."
	      sys.exit(1)

	# sprawdzam, czy sekondary strakczer jest przewidziana tak samo, jak sie rozni to dajemy [C]oil
	struktura_wynikowa = ""
	for j in range(len(structures[0])):
		wyraz = ""
		porownywacz= ""
		for i in range(len(structures)):
			wyraz = wyraz + structures[i][j]
			porownywacz = porownywacz + wyraz[0]
		if wyraz==porownywacz:
			struktura_wynikowa = struktura_wynikowa + wyraz[0]
		else:
			struktura_wynikowa = struktura_wynikowa + "C"

	print_seq_format(sequence, struktura_wynikowa)
