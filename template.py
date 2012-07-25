#!/usr/bin/env python
from re import compile,match
from numpy import reshape,linalg
class Template:
	
	"""
		Class used for template storage and atoms distance calculation
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
	def distance(self,idx_i,idx_j):
		i = int(idx_i)
		j = int(idx_j)
		if i in self.residues and j in self.residues:
			vi = self.coordinates[i]
			vj = self.coordinates[j]
			return linalg.norm(vi-vj)




a = Template("./playground/2pcy.pdb")
for i in range(1,98):
	print i,i+1, a.distance(i,i+1)
	
