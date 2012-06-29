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
		self.seqlen = len(sequence)
		
		self._createSEQ() # we always need SEQ file so I put it here


	def _createSEQ(self):
		"""
			Create SEQ input file for CABS, which contains sequence and (predicted) secondary structure of target chain.
			
			SEQ file is three-column file, where first column is index of residue (1...N), 
			second is three letter aminoacid, and third is secondary structure (1: coil, 2: helix, 4: beta)
			
			SEQ is required, and has to be created each time so we added it to __init__ of this class. 
			.. warning::
			It overrides any SEQ in working directory.
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
		for i in range(self.seqlen):
			if self.sequence[i] not in seqk:
				raise Errors("Problem with sequence letter "+self.sequence[i])
			if self.ss[i] not in ssk:
				raise Errors("Problem with secondary structure letter "+ self.ss[i]+ \
					" (only H - helix, E - extended/beta, C - coil/everything else")	
		# write SEQ file
		try:
			f = open("SEQ","w")
			frmt = "%5d   %3s    %1d\n"
			for i in range(self.seqlen):
				s = frmt % (i+1,seq_trans[self.sequence[i]],ss_trans[self.ss[i]])
				f.write(s)
			print Info("SEQ file created")
			f.close()
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
	def createLatticeReplicas(self,start_structures_fn=[],replicas=20):
		"""
			Create protein models projected onto CABS lattice, which will be used as replicas.
			
			:param start_structures_fn: list of paths to pdb files which should be used instead of templates models.  This parameter is optional, and probably not often used. Without it script creates replicas from templates files.
			:type start_structures_fn: list
			:param replicas: define number of replicas in CABS simulation. Usually you don't need to change it, default 20 is optimal for most cases. If number of replicas is smaller than number of templates - program will create replicas from first R templates. If there is less templates than replicas, replicas are generated sequentially from template models.
			:type replicas: integer
		"""
		if len(start_structures_fn)==0 and len(self.templates_fn)==0:
			raise Errors("lists start_structures_fn OR templates_filenames cannot be empty !")
			
		tempdir = mkdtemp('','CABStmp_','.')
		
		temp_filenames = self.templates_fn
		if len(start_structures_fn)>0:	# for user selected start models
			temp_filenames = start_structures_fn
			
		# create lattice model of all pdbs
	
		i = 0
		for tfn in temp_filenames:
			l = 0
			try:
				print Info("creating lattice model for "+tfn)
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
				arg = "../FF/a.out %d %d %d" % (self.seqlen,l,i) # TODO
				chainstart = Popen([arg], shell=True, stdout=PIPE)
				chainstart.communicate()
			
				chdir("../")
				i += 1
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				
		# create FCHAINS file which contains replicas coordinates
		
		m = len(temp_filenames) # number of templates
		chains_data = []
		for i in range(m): # load all chain files into memory
			chain = open(tempdir+"/CHAIN"+str(i),"r")
			chains_data.append(chain.readlines())
			chain.close()
			
	
		counter = 0
		k = 0 					# index of template
		fchains = open("FCHAINS","w")				
		while counter<replicas: 
			for line in chains_data[k]:
				fchains.write(line)
			k += 1
			if k==m: k=0
			counter +=1
		fchains.close()	
		print Info("FCHAINS file created")
			

				
			

class Info():
	"""
		Simple message system
	"""
	def __init__(self,text):
		self.text = text
	def __str__(self):
		return "INFO: "+self.text

		
class Errors(Exception):
	"""
		Simple error messages
	"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)



# tests
if __name__ == "__main__":
	seq = "IDVLLGADDGSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDASKISMSEEDLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTVN"
	ss =  "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEEEEEEEEEEEEEEEEEEEEEEE"
	templates = ["playground/2pcy_CA.pdb","playground/2pcy.pdb"]
	a = CABS(seq,ss,templates)

	a.createLatticeReplicas()
