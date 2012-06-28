#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tempfile import mkdtemp
from re import search
from os import getcwd,chdir
from subprocess import Popen,PIPE

class CABS:
	"""CABS main class.
	
	:param sequence: one line sequence of the target protein
	:type sequence: string
	:param secondary_structure: one line secondary structure for the target protein
	:type secondary_structure: string
	:param templates_filenames: path to 3D protein model templates in pdb file format which you want to use for modeling. Cα numbering in templates must be aligned to target sequence
	:type templates_filenames: list
	"""

	def __init__(self,sequence,secondary_structure,templates_filenames):
		if len(sequence)!=len(secondary_structure):
			raise Errors("Different lenght of sequence and secondary structure")
		self.sequence = sequence
		self.ss = secondary_structure
		self.templates_fn = templates_filenames

	def createSEQ(self):
		"""
			Create SEQ input file for CABS, which contains sequence and (predicted) secondary structure of target
		"""
		seq_trans={'A': "ALA",'R': "ARG",'N':"ASN",'D':"ASP",'C':"CYS",\
		'E':"GLU",'Q':"GLN",'G':"GLY",'H':"HIS",'I':"ILE",'L':"LEU",\
		'K':"LYS",'M':"MET",'F':"PHE",'P':"PRO",'S':"SER",'T':"THR",\
		'W':"TRP",'Y':"TYR",'V':"VAL"}
		ss_trans={'C':1,'H':2,'E':4} 
		# secondary structure 1=coil, 2=helix,  4=beta

		seqk = seq_trans.keys()
		ssk = ss_trans.keys()
		# check if seq/ss are correct
		for i in range(len(self.sequence)):
			if self.sequence[i] not in seqk:
				raise Errors("Problem with sequence letter "+self.sequence[i])
			if self.ss[i] not in ssk:
				raise Errors("Problem with secondary structure letter "+ self.ss[i]+ \
					" (only H - helix, E - extended/beta, C - coil/everything else")	
		# write SEQ file
		try:
			f = open("SEQ","w")
			frmt = "%5d   %3s    %1d\n"
			for i in range(len(self.sequence)):
				s = frmt % (i+1,seq_trans[self.sequence[i]],ss_trans[self.ss[i]])
				f.write(s)
			f.close()
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
	def createLatticeModel(self,start_structure_fn=''):
		"""
			create model projected onto CABS lattice
		"""
		tempdir = mkdtemp('','CABS_','.')

		if start_structure_fn != '':	# for user selected start model
			l=0 # template length
			try: 
				f = open(start_structure_fn)
				fw = open(tempdir+"/ALIGN0","w")
				
				for line in f.readlines():
					if search("^ATOM.........CA", line): # get only Calpha atoms
						fw.write(line)
						l +=1
				fw.close()		
				f.close()
				chdir(tempdir)
				
				arg = "../FF/a.out %d %d %d" % (len(self.sequence),l,0) # TODO
				chainstart = Popen([arg], shell=True, stdout=PIPE)
				chainstart.communicate()
				chdir("../")
			
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
		else:
			# create lattice model for all templates
			sqlen = len(self.sequence)
			
			i = 0
			for tfn in self.templates_fn:
				l = 0
				try:
					f = open(tfn,"r")
					output = "/ALIGN%d" % (i)

					fw = open(tempdir+output,"w")
					for line in f.readlines():
						if search("^ATOM.........CA", line): # get only Calpha atoms
							fw.write(line)
							l +=1
					fw.close()		
					f.close()
					
					chdir(tempdir)
					arg = "../FF/a.out %d %d %d" % (sqlen,l,i) # TODO
					chainstart = Popen([arg], shell=True, stdout=PIPE)
					chainstart.communicate()
					chdir("../")
					i += 1
				except IOError as e:
					print "I/O error({0}): {1}".format(e.errno, e.strerror)
				
				
				
			


class Errors(Exception):
	"""
		Simple error messages
	"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)



# tests

seq = "IDVLLGADDGSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDASKISMSEEDLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTVN"
ss =  "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEEEEEEEEEEE"
templates = ["playground/2pcy_CA.pdb","playground/2pcy.pdb"]
a = CABS(seq,ss,templates)
a.createSEQ()
a.createLatticeModel()
