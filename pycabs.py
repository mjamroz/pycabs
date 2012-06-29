#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tempfile import mkdtemp
from re import search,sub
from os import getcwd,chdir
from subprocess import Popen,PIPE
from shutil import copyfile
from time import clock

class CABS:
	"""CABS main class.
	
	:param sequence: one line sequence of the target protein
	:type sequence: string
	:param secondary_structure: one line secondary structure for the target protein
	:type secondary_structure: string
	:param templates_filenames: path to 3D protein model templates in pdb file format which you want to use for modeling. Cα numbering in templates must be aligned to target sequence
	:type templates_filenames: list
	:param project_name: files for that particular modeling will be named with prefix project_name.
	:type project_name: string
	
	.. warning:: it will overwrite previous SEQ file.
	"""

	def __init__(self,sequence,secondary_structure,templates_filenames,project_name="fnord"):
		if len(sequence)!=len(secondary_structure):
			raise Errors("Different lengths of sequence and secondary structure")
		self.sequence = sequence
		sub(r'\s', '', project_name)
		self.name = project_name
		self.ss = secondary_structure
		self.templates_fn = templates_filenames
		self.seqlen = len(sequence)
		
		self.constraints = 0
		# self.replicas - defined in lattie model creation method
		self._createSEQ() # we always need SEQ file so I put it here

	def calcConstraints(self,exclude_residues=[]): # TODO - update self.constraints
		pass
	def trafToPdb(self, output_filename = "TRAF.pdb"):
		pass	
	def _copyFFFiles(self):
		path_to_ff_dir="/home/hydek/pycabs/FF" # TODO
		for fff in ["R13","R13A","CENTRO","QUASI3S","R13C","R13E",\
		"R13H","R14","R14A","R14C","R14E","R14H","R15","R15A","R15C",\
		"R15E","R15H","SIDECENT"]:
			copyfile(path_to_ff_dir+"/"+fff,"./"+fff)
			
	def copyChains(self):
		copyfile("ACHAINS_NEW","FCHAINS")
		
	def modeling(self, Ltemp=1.0,Htemp=2.0,cycles=20,phot=10,constraints_force=1.0):
		self._copyFFFiles()
		
		# create INP file
		rng_seed = 1799
		inp = "%d %d %d %d\n%5.2f %5.2f %5.2f %5.3f %5.3f\n%5.2f %5.2f"\
		 "%5.2f %5.2f %5.3f\n0. 0 0 0 0\n0.0 0.0 0.0 0.0\n0 0.0\n%6d"\
		 "%5.2f" %(rng_seed,cycles,phot,self.replicas,Htemp,Ltemp,4.0,\
		 0.125,0.250,1.0,2.0,0.25,-2.00,0.375,self.constraints,constraints_force)
		try:
			f = open("INP","w")
			f.write(inp)
			f.close()
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
		print Info("INP file created")	
		# run CABS	
		# TODO zmodyfikowac CABS odnosnie OUT i TRAF
		print Info("CABS modeling software...")
		path_to_cabs="/home/hydek/pycabs/FF/" # TODO
		arg = path_to_cabs+"cabs" 
		cabs_start = clock()

		cabsstart = Popen([arg], shell=True, stdout=PIPE)
		cabsstart.communicate()

		elapsed = (clock() - cabs_start)/60

		print Info("Done. Running time: %6d min." %(elapsed))
		
		
	def _createSEQ(self):
		"""
			Create SEQ input file for CABS, which contains sequence and (predicted) secondary structure of target chain.
			
			SEQ file is three-column file, where first column is index of residue (1...N), 
			second is three letter aminoacid, and third is secondary structure (1: coil, 2: helix, 4: beta)
			
			SEQ is required, and has to be created each time so we added it to __init__ of CABS class. 
			**WARNING**: It overwrites any SEQ in working directory.
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
			:param replicas: define number of replicas in CABS simulation. However 20 is optimal for most cases, and you don't need to change it in protein modeling case. 
			:type replicas: integer
			
			.. note:: If number of replicas is smaller than number of templates - program will create replicas using first *replicas* templates. If there is less templates than replicas, they are creating sequentially using template models.
			
			.. warning:: it will overwrite FCHAINS file from working directory
		"""
		if len(start_structures_fn)==0 and len(self.templates_fn)==0:
			raise Errors("lists start_structures_fn OR templates_filenames cannot be empty !")
			
		tempdir = mkdtemp('',self.name+'_tmp','.')
		self.tempdir = tempdir
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
				arg = "/home/hydek/pycabs/FF/a.out %d %d %d" % (self.seqlen,l,i) # TODO
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
		self.replicas = replicas
		print Info("FCHAINS file created")
			

				
			


#33############ Utilities ########################
def parsePsipredOutput(psipred_output_fn):
	"""
	Psipred (protein secondary structure prediction, http://bioinf.cs.ucl.ac.uk/psipred/) output parser. 
	Psipred output looks like: ::
	
			> head psipred.ss
			1 P C   1.000  0.000  0.000
			2 K C   0.665  0.000  0.459
			3 A E   0.018  0.000  0.991
			4 L E   0.008  0.000  0.997
			5 I E   0.002  0.000  0.998
			6 V E   0.003  0.000  0.999
			7 Y E   0.033  0.000  0.981

	:param psipred_output_fn: path to the psipred output file
	:type psipred_output_fn: string
	:returns: tuple (sequence, secondary_structure)
		
	"""
	coil = ["-","T","S","B","G","I"]
	seq = ""
	ss = ""
	try:
		f = open(psipred_output_fn)
		for line in f.readlines():
			d = line.split()
			if d[2] in coil: d[2] = "C"
			seq = seq + d[1]
			ss = ss + d[2]
		f.close()
	except IOError as e:
		print "I/O error({0}): {1}".format(e.errno, e.strerror)
	return (seq,ss)
	
def parsePorterOutput(porter_output_fn):
	"""
		Porter (protein secondary stucture prediction, http://distill.ucd.ie/porter/) output parser. 
		Porter output looks like: ::
		
			IDVLLGADDGSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDASKISMSEE
			CEEEECCCCCCCCEECCEEEECCCCEEEEEECCCCCEEEEECCCCCCCCCCHHHHCCCCC



			DLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTVN
			CCECCCCCEEEEECCCCEEEEEECCHHHHCCCEEEEEEC

		.. note:: Use only sequence/secondary structure fragment of Porter email
		
		:param porter_output_fn: path to the porter output file
		:type porter_output_fn: string
		:returns: tuple (sequence, secondary_structure)
		
	"""
	coil = ["-","T","S","B","G","I"]
	seq = ""
	ss = ""
	try:
		f = open(porter_output_fn)
		i=0
		for line in f.readlines():
			line = line.strip()
			if line=="":
				continue
			if (i%2)==0:
				seq = seq+line
			else:
				ss = ss+line	
			i += 1
		f.close()
		# replace other than HCE into C
		for c in coil: ss = ss.replace(c,'C')
	except IOError as e:
		print "I/O error({0}): {1}".format(e.errno, e.strerror)
	return (seq,ss)
	



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


#################################################################




# tests
if __name__ == "__main__":
	data =  parsePorterOutput("playground/porter.ss")

	templates = ["playground/2pcy_CA.pdb"]
	a = CABS(data[0],data[1],templates)
	a.createLatticeReplicas()
	a.modeling()
#	print parsePsipredOutput("playground/psipred.ss")
