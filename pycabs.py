#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
pyCABS Copyright (C) 2012 Michal Jamroz <jamroz@chem.uw.edu.pl>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""


from tempfile import mkdtemp
from re import search,sub,compile,match
from os import getcwd,chdir,path,mkdir,stat,remove,rename
from numpy import fromfile,reshape,linalg,mean,std,zeros,unique,indices
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
		self.FF = "/home/mjamroz/pycabs/FF" # TODO !!!
		if len(sequence)!=len(secondary_structure):
			raise Errors("Different lengths of sequence and secondary structure")
		if path.isdir(project_name):
			raise Errors("Directory "+project_name+" exists. Choose another project name")
		self.sequence = sequence
		self.seqlen = len(sequence)
		
		self.ss = secondary_structure
		if len(templates_filenames)==0: # if 0 templates = DE NOVO 
			templates_filenames.append(self.__getInitForDeNovo())
			print templates_filenames
			
		self.templates_fn = templates_filenames
		self.cwd = getcwd()
		sub(r'\s', '', project_name)
		self.pname = project_name
		mkdir(self.pname)
		chdir(self.pname)
		
		self.rng_seed = 1799	#: seed for random generator
		self.constraints = 0
		
		# self.replicas - defined in lattice model creation method
		
		
		self.seq_trans={'A': "ALA",'R': "ARG",'N':"ASN",'D':"ASP",'C':"CYS",\
		'E':"GLU",'Q':"GLN",'G':"GLY",'H':"HIS",'I':"ILE",'L':"LEU",\
		'K':"LYS",'M':"MET",'F':"PHE",'P':"PRO",'S':"SER",'T':"THR",\
		'W':"TRP",'Y':"TYR",'V':"VAL"}
		
		self._createSEQ() # we always need SEQ file so I put it here
	def __getInitForDeNovo(self):
		
		f_chain =open(self.FF+"/extended.pdb").readlines()
		fw = open("start.pdb","w")
		fw.write("".join(f_chain[:self.seqlen]))
		fw.close()
		return path.abspath("start.pdb")
		
		
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

	def generateConstraints(self,exclude_residues=[],other_constraints=[]): # TODO - update self.constraints
		"""
			Calculate distance constraints using templates 3D models. Constraint will be a square well of size d-std_dev-1.0,d+std_dev+1.0, where d is mean distance among templates between Cα atoms (if constraint will be exceeded, there is penalty, scaled by weight.
			
			Weight is defined as a fraction of particular average distance among templates i.e. if pair of residues exist in 2 of 3 templates, weight will be 0.66. Using multiple sequence alignments it should provide stronger constraints on consistently aligned parts. 
			
			:param exclude_residues: indexes of residues without constrains
			:type exclude_residues: list
			:param other constrains: user-defined constrains as list of tuples: (residue_i_index,residue_j_index,distance, constraint_strength)
			:type other_constraints: list
			
		"""
		t = [] # here is list with templates structures
		for template_path in self.templates_fn:
			model = Template(template_path)
			t.append(model)
		max_resid = -1
		for temp in t:
			if temp.getLastResidueIndex()>max_resid:
				max_resid = temp.getLastResidueIndex()
		
		constraints = []
		for i in range(1,max_resid):
			later_j = []
			jtmp  = 14
			for l in range(14):
				jtmp = int(jtmp*1.4)
				if jtmp+i<max_resid:
					later_j.append(jtmp)
				else:
					break
			for j in [5,14]+later_j:
				ji = j+i
				if ji>max_resid or i in other_constraints or ji in other_constraints or i in exclude_residues or ji in exclude_residues:
					continue
					
				distances = []
				for temp in t:
					d = temp.distance(i,ji)
					if d:
						distances.append(d)
				if(len(distances)>0):
					d_mean = mean(distances)
					d_stddev = std(distances)
					d_weight = 1.0*len(distances)/len(t)
					
					d_min = d_mean-d_stddev-1.0
					if d_min<3.8: d_min = 3.8
					d_max = d_mean+d_stddev+1.0
					constraint = (i,ji,d_min,d_max,d_weight)
					constraints.append(constraint)
		self.constraints = len(constraints+other_constraints)
		f = open("constraints.dat","w")
		iter = 1
		for c in constraints+other_constraints:
			f.write( "%4d %5d %5d %6.2f %6.2f %5.2f\n"%(c[0],c[1],iter,round(c[2],2),round(c[3],2),round(c[4],2)))
			iter += 1
		f.close()
						
	def generateConstraintsOld(self,exclude_residues=[],other_constraints=[]):			
		"""
			Calculate distance constraints using templates 3D models. Constraint will be a square well of size min(d), max(d) where d is mean distance among templates between Cα atoms (if constraint will be exceeded, there is penalty, scaled by weight.
			
			Weight is defined as a fraction of particular average distance among templates i.e. if pair of residues exist in 2 of 3 templates, weight will be 0.66. Using multiple sequence alignments it should provide stronger constraints on consistently aligned parts. 
			
			:param exclude_residues: indexes of residues without constrains
			:type exclude_residues: list
			:param other constrains: user-defined constrains as list of tuples: (residue_i_index,residue_j_index,distance, constraint_strength)
			:type other_constraints: list
			
		"""
		t = [] # here is list with templates structures
		for template_path in self.templates_fn:
			model = Template(template_path)
			t.append(model)
		max_resid = -1
		for temp in t:
			if temp.getLastResidueIndex()>max_resid:
				max_resid = temp.getLastResidueIndex()
		
		constraints = []
		for i in range(1,max_resid):
			later_j = []
			jtmp  = 14
			for l in range(14):
				jtmp = int(jtmp*1.7)
				if jtmp+i<max_resid:
					later_j.append(jtmp)
				else:
					break
			for j in [5,14]+later_j:
				ji = j+i
				if ji>max_resid or i in other_constraints or ji in other_constraints or i in exclude_residues or ji in exclude_residues:
					continue
					
				distances = []
				for temp in t:
					d = temp.distance(i,ji)
					if d:
						distances.append(d)
				if(len(distances)>0):
					
					d_min = min(distances)
					d_max = max(distances)
					d_weight = len(distances)/(1.0+abs(d_max-d_min))
					
					constraint = (i,ji,d_min,d_max,d_weight)
					constraints.append(constraint)
		self.constraints = len(constraints+other_constraints)
		f = open("constraints.dat","w")
		iter = 1
		for c in constraints+other_constraints:
			f.write( "%4d %5d %5d %6.2f %6.2f %5.2f\n"%(c[0],c[1],iter,round(c[2],2),round(c[3],2),round(c[4],2)))
			iter += 1
		f.close()
        			
	def getEnergy(self):
		"""
			Read CABS energy values into list
			
			:return: list of models energy
			
		"""
		if path.isfile("ENERGY"):
			try:
				energy = fromfile(e_path,sep='\n') # read ENERGY data into array `energy`
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no ENERGY file in current directory, did you run CABS.modeling method before?")
			return energy
			
	def savePdbModel(self, model_idx, filename=''):
		"""
			Save trajectory model into pdb file
			
			:param model_idx: index of model in CABS trajectory
			:param filename: name of the output file. If empty, it saves to model_index.pdb
		"""
		if path.isfile("TRAF"):
			seq = []
			for s1 in self.sequence: seq.append(self.seq_trans[s1])
				
			try:
				traf = open("TRAF")
				t = traf.readlines()
				traf.close()
				model = []
				index = -1
				for l in t:
					if "." in l:
						index +=1
					if index==model_idx:
						for coord in l.split(): model.append(int(coord)*0.61)
					
				model = model[3:-3]
				pdb_format = "ATOM   %4d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
				model_format = "MODEL%9d\n"
				pdb_data = model_format %(model_idx)
				for ai in range(self.seqlen):
					pdb_data += pdb_format %(ai+1,seq[ai],ai+1,model[3*ai],model[3*ai+1],model[3*ai+2])
				pdb_data += "ENDMDL\n"		
				if filename=='':
					out = "model_%04d.pdb" % (model_idx)
				else:
					out = filename
					
				fw = open(out,"w")
				fw.write(pdb_data)
				fw.close()
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no TRAF file in current directory, did you run CABS.modeling method before?")				
			
			
	def getTraCoordinates(self):
		"""
			Read trajectory file into 2D list of coordinates
			
			:return: 2D list of trajectory coordinates ( list[1][6] is sixth coordinate of second trajectory model = z coordinate of second atom of second model)
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
				
	def convertPdbToDcd(self,catdcd_path="/home/user/pycabs/FF/catdcd"):
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
		
		
	def modeling(self, Ltemp=1.0,Htemp=2.0,cycles=1,phot=1,constraints_force=1.0,dynamics=False):
		"""
			Start CABS modeling
			
			:param Ltemp: Low temperature for Replica Exchange Monte Carlo
			:type Ltemp: float
			:param Htemp: High temperature for Replica Exchange Monte Carlo
			:type Htemp: float
			:param cycles: number of Replica Exchange cycles
			:type cycles: integer
			:param iphot: number of microcycles (inside REMC loop)
			:type iphot: integer
			:param constraints_force: Slope of constraints force potential
			:type constraints_force: float
            :param dynamics: Use of special CABS version for dynamics pathway studies
			:type dynamics: boolean
            
			
		"""
		#preprocessing
		self._copyFFFiles()
		
		# create INP file
		#TODO
		inp = "%5d %5d %5d %5d\n%5.2f %5.2f %5.2f %5.3f %5.3f\n%5.2f %5.2f"\
		 "%5.2f %5.2f %5.3f\n0. 0 0 0 0\n0.0 0.0 0.0 0.0\n0 0.0\n%5d"\
		 "%5.2f\n" %(self.rng_seed,cycles,phot,self.replicas,Htemp,Ltemp,4.0,\
		 0.125,0.250,1.0,2.0,0.25,-2.00,0.375,self.constraints,constraints_force)
		try:
			f = open("INP","w")
			f.write(inp)
			if(self.constraints>0):
				fc = open("constraints.dat")
				for line in fc.readlines():
					f.write(line)
				fc.close()
			
			f.close()
		except IOError as e:
			print "I/O error({0}): {1}".format(e.errno, e.strerror)
		print Info("INP file created")	
		
		# run CABS	
		print Info("CABS started...")
		arg = path.join(self.FF,"cabs")
        	if dynamics:
	            arg = path.join(self.FF,"cabs_dynamics")
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
	def loadSGCoordinates(self):
		"""
			Read center of mass of sidegroups from TRASG file

			:return: 2D list of sidechains coordinates

		"""
		filename = "TRASG"
		trajectory = []
		if path.isfile(filename):
			try:
				traf = open(filename)
				t = traf.readlines()
				traf.close()
				chain_len = int(t[0].split()[1])+1
				for i in range(0,len(t),chain_len):
					model = []
					for j in range(i+1,i+chain_len):
						coo = t[j].split()
						model.append(float(coo[1])*0.61)
						model.append(float(coo[2])*0.61)
						model.append(float(coo[3])*0.61)
					trajectory.append(model)
					
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no TRASG file in current directory, did you run CABS.modeling method before?")		
		return trajectory		
			
			



#33############ Utilities ########################
def parseDSSPOutput(filename):
    """
        Helper function for extracting sequence and secondary structure assignments from the DSSP output. Useful for dynamics studies or other where we know protein structure.
        
        You can download DSSP files directly from PDB server: http://www.pdb.org/pdb/files/PDBID.dssp
    """
    h = open(filename,"r")
    # translate cystein
    from string import maketrans
    tr = maketrans("abcdefghijklmnopqrstuvwxyz","CCCCCCCCCCCCCCCCCCCCCCCCCC")
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
		
    return (seq.translate(tr),ss)
    
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

def parsePDBfile(pdb_filename):
	"""
	
		Function for parsing of Cα coordinates from PDB file. 
		
		:param pdb_filename: path to PDB file
		:type pdb_filename: string
		:return: 1D list of Cα coordinates (for example: list[4] is y-th coordinate of second atom)
	
	"""
	f = open(pdb_filename)
	atm = compile(r"^ATOM.{9}CA.{7}(?P<resid>.{4})(?P<x>.{12})(?P<y>.{8})(?P<z>.{8})")
	model = []
	for line in f.readlines():
		data = atm.match(line)
		if data:
			for v in data.groups()[1:]: model.append(float(v))
	f.close()
	return model
			
		
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
def heat_map(data, x_label,y_label,colormap_label,output_file="heatmap.png",cmap='Greys'):
	"""
		Save heat map using pylab
		
		:param data: 2D list of values
		:type data: float
		
	"""
	
	from pylab import xlabel,ylabel,pcolor,colorbar,savefig 
	l = len(data[0])
	rows, cols = indices((l,l))
	xlabel(x_label)
	ylabel(y_label)
	pcolor(data, cmap=cmap)
	cb = colorbar()
	cb.set_label(colormap_label)
	savefig(output_file)

def loadTRAFCoordinates(filename):
		"""
			Read trajectory file into 2D list of coordinates
			
			:return: 2D list of trajectory coordinates ( list[1][6] is sixth coordinate of second trajectory model = z coordinate of second atom of second model)
		"""

		trajectory = []
		if path.isfile(filename):
			try:
				traf = open(filename)
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

def loadSGCoordinates(filename):
		"""
			Read center of mass of sidegroups from TRASG file

			:param filename: path to the TRASG file
			:return: 2D list of sidechains coordinates

		"""
		trajectory = []
		if path.isfile(filename):
			try:
				traf = open(filename)
				t = traf.readlines()
				traf.close()
				chain_len = int(t[0].split()[1])+1
				for i in range(0,len(t),chain_len):
					model = []
					for j in range(i+1,i+chain_len):
						coo = t[j].split()
						model.append(float(coo[1])*0.61)
						model.append(float(coo[2])*0.61)
						model.append(float(coo[3])*0.61)
					trajectory.append(model)
					
			except IOError as e:
				print "I/O error({0}): {1}".format(e.errno, e.strerror)
				raise Errors("Maybe there is no TRASG file in current directory, did you run CABS.modeling method before?")		
		return trajectory
		
def saveMedoids(clusters,cabs):
	"""
		Save cluster medoids in PDB file format.
		
		:param clusters: cluster indices as a output of C Clustering Library
		:param clusters: list
		
	"""
	seq = []
	for s1 in cabs.sequence: seq.append(cabs.seq_trans[s1])
	pdb_format = "ATOM   %4d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n"
	trajectory = cabs.getTraCoordinates()
	
	medoids_idx = []
	for c in unique(clusters):
		centroid = zeros((len(trajectory[0])))
		this_cluster_idx = []
		cluster_size = 0
		for i in range(len(clusters)):
			if clusters[i]==c:
				cluster_size +=1
				model = array(trajectory[i])
				this_cluster_idx.append(i)
				centroid = centroid + model
		centroid=centroid/cluster_size
		min_rms = 50000
		min_idx = -1
		for e in this_cluster_idx:
			rms = rmsd(centroid,trajectory[e])
			if rms<min_rms:
				min_rms = rms
				min_idx = e
		medoids_idx.append(min_idx)
	for i in range(len(medoids_idx)):
		pdb_data = "MODEL%9d\n" %(medoids_idx[i])
		outfile = "model_%04d.pdb" %(i+1)
		model = trajectory[medoids_idx[i]]
		for ai in range(cabs.seqlen):
			pdb_data += pdb_format %(ai+1,seq[ai],ai+1,model[3*ai],model[3*ai+1],model[3*ai+2])
		pdb_data += "ENDMDL\n"
		fw = open(outfile,"w")
		fw.write(pdb_data)
		fw.close()
					
	
	return
def contact_map(trajectory, contact_cutoff):
	"""
		Compute fraction of contacts in trajectory, where trajectory is 2D list of coordinates (trajectory[2][5] is the z-th coordinate of second atom of third model)

		:param trajectory: 2D trajectory of atoms (Cα, sidegroups center of mass, etc.)
		:param contact_cutoff: cutoff defining contact

		:return: 2D array of fraction of contacts (number of contacts divided by trajectory length) for each pair of residue.
	"""

	model_len = len(trajectory[0])/3 # model is a 1D list of coordinates
	contacts_tmp = zeros((model_len,model_len))
	cutoff2=contact_cutoff*contact_cutoff
	for model in trajectory:
		for i in range(model_len-1):
			for j in range(i+1,model_len):
				x_dist = model[3*i] - model[3*j]
				y_dist = model[3*i+1] - model[3*j+1]
				z_dist = model[3*i+2] - model[3*j+2]
				d = sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)
				if d<contact_cutoff:
					contacts_tmp[i][j] += 1.0
					contacts_tmp[j][i] += 1.0
	#contacts_tmp/=len(trajectory)
	return contacts_tmp			
				
	

		
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
            
			:return: array of model coordinates
            
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
   
		:param reference: 1D list of coordinates (length of 3N)
		:type reference: list
		:param arr: 1D list of coordinates (length of 3N)
		:type arr: list
		:return: RMSD value after optimal superimposition of two structures
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
    
	if (root0<1e-6): root0 = 0.0
	if (root1<1e-6): root1 = 0.0
	if (root2<1e-6): root2 = 0.0
    
	roots = [root0,root1,root2]

	roots.sort()
	minr=roots[0]
	midr=roots[1]
	maxr=roots[2]
	dwa = 2.0*(sqrt(maxr) + sqrt(midr) + sign*sqrt(minr) )
	if dwa>Rg:
		rms = 0.0
	else:
		rms = sqrt(Rg - dwa)
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

class Template:
	
	"""
		Class used for storage of templates atom positions and distance calculation
		
		:param filename: path to file with template (in PDB format)
		:return: Nx3 list of coordinates
	"""

	def __init__(self,filename):
		
		self.hashmap={}
		atm = compile(r"^ATOM.{9}CA.{7}(?P<resid>.{4})(?P<x>.{12})(?P<y>.{8})(?P<z>.{8})")
		f = open(filename)
		coordinateslist = []
		i=0
		for line in f.readlines():
			data = atm.match(line) # parse pdb file
			if data:
				self.hashmap[int(data.group('resid'))] = i # remember residues indexes
				i+=1
				for v in data.groups()[1:]: coordinateslist.append(float(v))
		f.close()
		self.coordinates = reshape(coordinateslist,(-1,3)) # magic, reshape list of xyzxyz into Nx3 array
		coordinateslist = []
		self.residues = self.hashmap.keys()
	def printCoordinates(self):	
		for i in self.hashmap.keys():
			print i,self.hashmap[i],self.coordinates[self.hashmap[i]]
	def getCoordinate(self,residue_index):
		idx = int(residue_index)
		if idx in self.residues:
			return self.coordinates[idx]
	def getLastResidueIndex(self):
		return self.residues[-1]
	def distance(self,idx_i,idx_j):
		"""
			:param idx_i: residue index (as in target sequence numbering)
			:param idx_j: residue index (as in target sequence numbering)
			:type idx_i: integer
			:type idx_j: integer
			:return: euclidean distance between Cα(i) and Cα(j)
		"""
		i = int(idx_i)
		j = int(idx_j)
		if i in self.residues and j in self.residues:
			vi = self.coordinates[self.hashmap[i]]
			vj = self.coordinates[self.hashmap[j]]
			return linalg.norm(vi-vj)
		else:
			return None

		
#################################################################




# tests
if __name__ == "__main__":
	#data =  parsePorterOutput("/home/user/pycabs/proba/playground/porter.ss") # read PORTER (or PsiPred) secondary structure prediction
	seq = "IDVLLGADDGSLAFVPSEFSISPGEKIVFKNNAGFPHNIVFDEDSIPSGVDASKISMSEEDLLNAKGETFEVALSNKGEYSFYCSPHQGAGMVGKVTVN"
	ss = "CEEEECCCCCCCEEECCEEEECCCCEEEEEECCCCCCEEEECCCCCCCCCCCCCCCCCCCCCCCCCCCEEEEECCCCEEEEEECCCCCCCCCEEEEEEC"
	working_dir = "modelowanie2pcy" # name of project 
	templates = ["/home/mjamroz/pycabs/playground/2pcy_CA.pdb"]#,"/home/mjamroz/pycabs/playground/2pcy_CA2.pdb"] # set path to templates 
	a = CABS(seq,ss,templates,working_dir) # initialize CABS, create required files
	a.generateConstraintsOld() 
	#a.createLatticeReplicas() # create start models from templates
	#a.modeling(Htemp=3.0,cycles=1,phot=1) # start modeling with default INP values and create TRAF.pdb when done
	#tr = a.getTraCoordinates() # load TRAF into memory and calculate RMSD all-vs-all : 
	
	clu = '''
	from Pycluster import *
	from numpy import array
	distances = zeros((len(tr),len(tr)))
	for i in range(len(tr)):
		for j in range(i,len(tr)):
			rms = rmsd(tr[i],tr[j])
			distances[i][j] = distances[j][i] = rms
	heat_map(distances,"Protein model","Protein model","RMSD")		
	
	clusterid,error,nfound = kmedoids(distances,nclusters=5,npass=15,initialid=None)
	print clusterid,error
	clusterid,error,nfound = kcluster(distances,nclusters=5,npass=15)
	saveMedoids(clusterid,a)
	print clusterid,error
    	'''
