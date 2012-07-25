#!/usr/bin/env python
from re import compile,match
from numpy import reshape
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
	def printCoordinates(self):	
		for i in self.hashmap.keys():
			print i,self.hashmap[i],self.coordinates[self.hashmap[i]]





a = Template("./playground/2pcy.pdb")
a.printCoordinates()
