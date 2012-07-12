#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tempfile import mkdtemp
from re import search,sub
from os import getcwd,chdir,path,mkdir,stat,remove,rename
from subprocess import Popen,PIPE
from shutil import copyfile
from math import sqrt,cos,sin,atan2
import threading
import time

class CABS(threading.Thread):
	"""
	
	CABS main class.
	
	:param sequence: one line sequence of the target protein
	:type sequence: string
	:param secondary_structure: one line secondary structure for the target protein
	:type secondary_structure: string
	:param templates_filenames: path to 3D protein model templates in pdb file format which you want to use for modeling. Cα numbering in templates must be aligned to target sequence
	:type templates_filenames: list
	:param project_name: project_name and working directory name (uniq)
	:type project_name: string
	
	"""

	def __init__(self,sequence,secondary_structure,templates_filenames, project_name):
		self.FF = "/home/hydek/pycabs/FF" # TODO !!!
		if len(sequence)!=len(secondary_structure):
			raise Errors("Different lengths of sequence and secondary structure")
		if path.isdir(project_name):
			raise Errors("Directory "+project_name+" exists. Choose another project name")
		self.sequence = sequence
		self.ss = secondary_structure
		self.templates_fn = templates_filenames
		self.cwd = getcwd()
		sub(r'\s', '', project_name)
		self.pname = project_name
		mkdir(self.pname)
		chdir(self.pname)
		
		self.seqlen = len(sequence)
		self.rng_seed = 1799	#: seed for random generator
		self.constraints = 0
		# self.replicas - defined in lattice model creation method
		
		
		self.seq_trans={'A': "ALA",'R': "ARG",'N':"ASN",'D':"ASP",'C':"CYS",\
		'E':"GLU",'Q':"GLN",'G':"GLY",'H':"HIS",'I':"ILE",'L':"LEU",\
		'K':"LYS",'M':"MET",'F':"PHE",'P':"PRO",'S':"SER",'T':"THR",\
		'W':"TRP",'Y':"TYR",'V':"VAL"}
		
		self._createSEQ() # we always need SEQ file so I put it here

	def calcConstraints(self,exclude_residues=[],other_constraints=[]): # TODO - update self.constraints
		"""
			Calculate distance constraints using templates 3D models. 
			
			:param exclude_residues: indexes of residues without constrains
			:type exclude_residues: list
			:param other constrains: user-defined constrains as list of tuples: (residue_i_index,residue_j_index,constraint_strength)
			:type other_constraints: list
			
		"""
		pass
	def getTraCoordinates(self):
		"""
			read trajectory file into 2D list of coordinates
			
		"""
		trajectory = []
		if path.isfile("TRAF"):
			try:
				traf = open("TRAF")
				t = traf.readlines()
				traf.close()
				model = []
				for l in t[1:]:
					if "." in l: # "." is only in TRAF header	
						model = model[3:-3] # skip dummy atoms
						trajectory.append(model)
						model = []
					else:	# load coordinates
						for coord in l.split(): model.append(int(coord)*0.61) 
				model = model[3:-3] # skip dummy atoms
				trajectory.append(model)
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no TRAF file in current directory, did you run CABS.modeling method before?")		
		return trajectory
		
	def trafToPdb(self, output_filename = "TRAF.pdb"):
		""" 
			Convert TRAF CABS pseudotrajectory file format into multimodel pdb
		"""
		
		pdb_format = "ATOM   %4d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
		model_format = "MODEL%9d\n"

		if path.isfile("TRAF"):
			try:
				traf = open("TRAF")
				t = traf.readlines()
				traf.close()
				
				# it is worth to create 3-digit list of residues
				seq = []
				for s1 in self.sequence: seq.append(self.seq_trans[s1])
				pdb_data = ""
				model = []
				model_idx = 1
				for l in t[1:]:
					if "." in l: # "." is only in TRAF header
						pdb_data += model_format %(model_idx)
						model = model[3:-3] # skip dummy atoms
						for ai in range(self.seqlen):
							pdb_data += pdb_format %(ai+1,seq[ai],ai+1,model[3*ai],model[3*ai+1],model[3*ai+2])
						pdb_data += "ENDMDL\n"	
						model = []
						model_idx +=1 
					else:	# load coordinates
						for coord in l.split(): model.append(int(coord)*0.61) 
				# and write last model and clean a little
				pdb_data += model_format %(model_idx)
				model = model[3:-3] # skip dummy atoms
				for ai in range(self.seqlen):
					pdb_data += pdb_format %(ai+1,seq[ai],ai+1,model[3*ai],model[3*ai+1],model[3*ai+2])
				pdb_data += "ENDMDL\n"		
				t = []
				fw = open(output_filename,"w")
				fw.write(pdb_data)
				fw.close()
				pdb_data = ""
				
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no TRAF file in current directory, did you run CABS.modeling method before?")
				
	def convertPdbToDcd(self,catdcd_path="/home/hydek/pycabs/FF/catdcd"):
		"""
			This is only simple wrapper to CatDCD software (http://www.ks.uiuc.edu/Development/MDTools/catdcd/), 
			could be usable since \*.dcd binary format is few times lighter than pdb, and many python libraries 
			(ProDy, MDAnalysis) use \*.dcd as trajectory input format.
			Before use, download CatDCD from http://www.ks.uiuc.edu/Development/MDTools/catdcd/ and modify catdcd_path.
			
		"""
		output_dcd="TRAF.dcd"
		input_pdb="TRAF.pdb"
		if not path.isfile(catdcd_path):
			raise Errors("You do not have catdcd. Download it from http://www.ks.uiuc.edu/Development/MDTools/catdcd/ , uncompress and give correct path. Buenos dias.")
		if path.isfile(input_pdb):
			print Info("Converting PDB into DCD...")
			catdcd = Popen([catdcd_path+" -o "+output_dcd+" -otype dcd -pdb "+input_pdb], shell=True, stdout=PIPE)
			catdcd.communicate()
			
		else:
			raise Errors("input pdb file does not exist!")
	def _checkOverlaps(self):
		f = open("OUT")
		for line in f.readlines():
			if search("overlaps",line) and int(line.split()[2])>0:
				print Info(" **WARNING** "+line.strip())
		f.close()
					
	def _copyFFFiles(self):

		for fff in ["R13","R13A","CENTRO","QUASI3S","R13C","R13E",\
		"R13H","R14","R14A","R14C","R14E","R14H","R15","R15A","R15C",\
		"R15E","R15H","SIDECENT"]:
			copyfile(path.join(self.FF,fff),fff)
	def _removeFFFiles(self):
		for fff in ["R13","R13A","CENTRO","QUASI3S","R13C","R13E",\
		"R13H","R14","R14A","R14C","R14E","R14H","R15","R15A","R15C",\
		"R15E","R15H","SIDECENT"]:
			remove(fff)
				
	def _copyChains(self):
		copyfile("FCHAINS","FCHAINS_old")
		rename("ACHAINS_NEW","FCHAINS")
	def setParameters(self,	 Ltemp=1.0,Htemp=2.0,cycles=20,phot=10,constraints_force=1.0):
		pass
	def modeling(self, Ltemp=1.0,Htemp=2.0,cycles=1,phot=1,constraints_force=1.0):
		
		#preprocessing
		self._copyFFFiles()
		
		# create INP file
		#TODO
		inp = "%5d %5d %5d %5d\n%5.2f %5.2f %5.2f %5.3f %5.3f\n%5.2f %5.2f"\
		 "%5.2f %5.2f %5.3f\n0. 0 0 0 0\n0.0 0.0 0.0 0.0\n0 0.0\n%5d"\
		 "%5.2f" %(self.rng_seed,cycles,phot,self.replicas,Htemp,Ltemp,4.0,\
		 0.125,0.250,1.0,2.0,0.25,-2.00,0.375,self.constraints,constraints_force)
		try:
			f = open("INP","w")
			f.write(inp)
			f.close()
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
		print Info("INP file created")	
		
		# run CABS	
		print Info("CABS started...")
		arg = path.join(self.FF,"cabs")
		cabsstart = Popen([arg], shell=True, stdout=PIPE)
		cabsstart.communicate()
		print Info("CABS Done.")
		
		# postprocessing
		self._checkOverlaps()
		self._copyChains()
		self._removeFFFiles()
		self.trafToPdb()
		
		
	def _createSEQ(self):
		"""
			Create SEQ input file for CABS, which contains sequence and (predicted) secondary structure of target chain.
			
			SEQ file is three-column file, where first column is index of residue (1...N), 
			second is three letter aminoacid, and third is secondary structure (1: coil, 2: helix, 4: beta)
			
			SEQ is required, and has to be created each time so we added it to __init__ of CABS class. 
			**WARNING**: It overwrites any SEQ in working directory.
		"""
		ss_trans={'C':1,'H':2,'E':4} # secondary structure 1=coil, 2=helix,  4=beta

		seqk = self.seq_trans.keys()
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
				s = frmt % (i+1,self.seq_trans[self.sequence[i]],ss_trans[self.ss[i]])
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
			
		"""
		
		if len(start_structures_fn)==0 and len(self.templates_fn)==0:
			raise Errors("lists start_structures_fn OR templates_filenames cannot be empty !")
		
		tempdir = mkdtemp('','tmp','')
		self.tempdir = tempdir
		
		temp_filenames = self.templates_fn
		if len(start_structures_fn)>0:	# for user selected start models
			temp_filenames = start_structures_fn
			
		# create lattice model of each pdb
	
		i = 0
		for tfn in temp_filenames:
			l = 0
			try:
				print Info("creating lattice model for "+tfn)
				f = open(tfn,"r")
				fw = open(path.join(self.cwd,self.pname,tempdir,"ALIGN%d" % (i)),"w")
				for line in f.readlines():
					if search("^ATOM.........CA", line): # get only Calpha atoms
						fw.write(line)
						l +=1
				fw.close()		
				f.close()
				chdir(path.join(self.cwd,self.pname,tempdir))
				arg = path.join(self.FF,"a.out")+ " %d %d %d" % (self.seqlen,l,i)
				chainstart = Popen([arg], shell=True, stdout=PIPE)
				chainstart.communicate()
				
				chdir(path.join(self.cwd,self.pname))
				i += 1
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				
		# create FCHAINS file which contains replicas coordinates
		
		m = len(temp_filenames) # number of templates
		chains_data = []
		for i in range(m): # load all chain files into memory
			chain = open(path.join(self.cwd,self.pname,tempdir,"CHAIN"+str(i)),"r")
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
		Porter emailed output looks like: ::
		
			IDVLLGADDGSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDASKISMSEE
			CEEEECCCCCCCCEECCEEEECCCCEEEEEECCCCCEEEEECCCCCCCCCCHHHHCCCCC



			DLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTVN
			CCECCCCCEEEEECCCCEEEEEECCHHHHCCCEEEEEEC
		
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
			if search("[a-z]",line): # ignore lines with small letters - headers/comments
				continue
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


		
class Monitor(threading.Thread):
	"""
		Class for monitoring of CABS output data. You can run it and dynamically update output arrays with calculated results.
		
		:param calculate: what to do with gathered data ? 
		:type calculate: :class:`Calculate`
	"""
		
	def __init__(self,filename,calculate):
		threading.Thread.__init__(self)
		self.daemon = False #: if True, it will terminate when script terminates		
		self.f = filename
		self.first_mtime = stat(filename).st_mtime
		self.kill = False
		self.fb = open(filename)
        
		self.buf = [i.strip() for i in self.fb] # get whole file
		self.position = self.fb.tell() # and save position of eof
		self.calc = calculate
		
		self.calc.calculate(self.buf) #update calculate object
		
	def terminate(self):
		"""
			Terminate monitor
		"""
		self.kill=True
		
	def run(self):
		"""
			Run monitor in background
		"""
		
		self.buf=[]
		while 1:
			if self.kill:
				return
			time.sleep(1)
			mtime = stat(self.f).st_mtime
			if mtime!=self.first_mtime:
				self.first_mtime = mtime
				self.fb.seek(self.position)
				for l in self.fb:
					self.buf.append(l.strip())
				self.position = self.fb.tell() #update EOF
				
				self.calc.calculate(self.buf) #update calculate object after each fortran FLUSH 
				self.buf=[]
					
class Calculate:
	"""
		Inherit if you want to process data used with :class:`Monitor` class.
		
		:param output: output array with calculated results
		:type output: array/list
	"""
	def __init__(self,output):
		self.out = output
	def calculate(self,data):
		self.out.append('_'.join(data))	
    
	def processTrajectory(self,data):
		"""
			Use it in `calculate` method if you parsing TRAF file, and want to calculate something on structure
            
			:result: array of model coordinates
            
		"""
    
		i=0
		model=[]
		mtmp = []
		for line in data:
			if "." in line :
				if(len(mtmp)>0):
					model.append(mtmp)
					mtmp = []
			else:
				for e in line.split():
					mtmp.append(0.61*int(e))
		if len(mtmp)>0:
			model.append(mtmp[3:-3])
		return model


def rmsd(reference,arr):
	"""
		Calculate coordinate Root Mean Square Deviation between two sets of coordinates.
		
		.. math::
		
			cRMSD = \sqrt{ \sum_{i=1}^N \|x_{i} - y_{i}\|^2 \over N}
   
		:param reference: 1D list of coordinates (length 3N)
		:type reference: list
		:param arr: 1D list of coordinates (length 3N)
		:type arr: list
		:return: RMSD after optimal superimposition
	"""
		
	l = len(arr)
	invlen = 3.0/l
	# move to 0,0,0
	x_cm=y_cm=z_cm=0.0
	ax_cm=ay_cm=az_cm=0.0
	
	for a in range(0,l-2,3):
		x_cm +=reference[a]
		y_cm +=reference[a+1]
		z_cm +=reference[a+2]
		ax_cm +=arr[a]
		ay_cm +=arr[a+1]
		az_cm +=arr[a+2]
	x_cm *=invlen
	y_cm *=invlen
	z_cm *=invlen
	ax_cm *=invlen
	ay_cm *=invlen
	az_cm *=invlen
	for a in range(0,l-2,3):
		reference[a] -=x_cm
		reference[a+1] -=y_cm
		reference[a+2] -=z_cm
		arr[a] -= ax_cm
		arr[a+1] -= ay_cm
		arr[a+2] -= az_cm
		
	covmat0=covmat1=covmat2=covmat3=covmat4=covmat5=covmat6=covmat7=covmat8=Rg=0
	for a in range(0,l-2,3):
		s_i_x = reference[a]
		s_i_y = reference[a+1]
		s_i_z = reference[a+2]
		s_j_x = arr[a]
		s_j_y = arr[a+1]
		s_j_z = arr[a+2]
		covmat0 += s_i_x*s_j_x
		covmat1 += s_i_y*s_j_x
		covmat2 += s_i_z*s_j_x
		covmat3 += s_i_x*s_j_y
		covmat4 += s_i_y*s_j_y
		covmat5 += s_i_z*s_j_y
		covmat6 += s_i_x*s_j_z
		covmat7 += s_i_y*s_j_z
		covmat8 += s_i_z*s_j_z
		Rg +=  s_i_x*s_i_x + s_j_x*s_j_x
		Rg +=  s_i_y*s_i_y + s_j_y*s_j_y
		Rg +=  s_i_z*s_i_z + s_j_z*s_j_z
	
	Rg *= invlen
	covmat0 *= invlen
	covmat1 *= invlen
	covmat2 *= invlen
	covmat3 *= invlen
	covmat4 *= invlen
	covmat5 *= invlen
	covmat6 *= invlen
	covmat7 *= invlen
	covmat8 *= invlen
	
	determinant = covmat0*(covmat4*covmat8 - covmat5*covmat7) - \
	covmat1*(covmat3*covmat8-covmat5*covmat6) + covmat2*\
	(covmat3*covmat7-covmat4*covmat6)
	sign = 1.0
	if (determinant<0.0): sign = -1.0
	 
	r0 = covmat0*covmat0 +    covmat3*covmat3 + covmat6*covmat6
	r1 = covmat0*covmat1 +    covmat3*covmat4 + covmat6*covmat7
	r2 = covmat0*covmat2 +    covmat3*covmat5 + covmat6*covmat8
	r4 = covmat1*covmat1 +    covmat4*covmat4 + covmat7*covmat7
	r5 = covmat1*covmat2 +    covmat4*covmat5 + covmat7*covmat8
	r8 = covmat2*covmat2 +    covmat5*covmat5 + covmat8*covmat8
	
	inv3 = 1.0/3.0
	root3 = 1.732050807568877294
	c0 = r0*r4*r8 + 2.0*r1*r2*r5 - r0*r5*r5 - r4*r2*r2 - r8*r1*r1
	c1 = r0*r4 - r1*r1 + r0*r8 - r2*r2 + r4*r8 - r5*r5
	c2 = r0+r4+r8
	c2Div3 = c2*inv3
	aDiv3 = (c1-c2*c2Div3)*inv3
	if(aDiv3>0.0): aDiv3 = 0.0
	mbDiv2 = 0.5*(c0+c2Div3*(2.0*c2Div3*c2Div3 - c1))
	q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3
	if (q>0.0): q=0.0

	magnitude = sqrt(-aDiv3)
	angle = atan2(sqrt(-1.0*q),mbDiv2)/3.0;
	sn=sin(angle)
	cs=cos(angle)
	magnitudecs = magnitude*cs
	magnituderoot = magnitude*root3*sn

	root0 = c2Div3 + 2.0*magnitudecs
	root1 = c2Div3 - magnitudecs - magnituderoot
	root2 = c2Div3 - magnitudecs +magnituderoot
    
	if (abs(root0)<1e-6): root0 = 0.0
	if (abs(root1)<1e-6): root1 = 0.0
	if (abs(root2)<1e-6): root2 = 0.0
    
	roots = [root0,root1,root2]

	roots.sort()
	minr=roots[0]
	midr=roots[1]
	maxr=roots[2]
	dwa = 2.0*(sqrt(maxr) + sqrt(midr) + sign*sqrt(minr) )
	rms = sqrt(Rg - dwa )
	if(rms<1e-5): rms=0.0
	return rms
	
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
	data =  parsePorterOutput("/home/hydek/pycabs/proba/playground/porter.ss")

	working_dir = "modelowanie2pcy"
	templates = ["/home/hydek/pycabs/playground/2pcy_CA.pdb"]
	a = CABS(data[0],data[1],templates,working_dir)
	a.createLatticeReplicas()
	a.modeling()
	tr = a.getTraCoordinates()
	print "RMSD ",rmsd(tr[1],tr[2])
	print tr[-1]
	a.convertPdbToDcd()
    
#	print parsePsipredOutput("playground/psipred.ss")
