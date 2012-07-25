#!/usr/bin/env python
from re import compile,match

class Template:
	"""
		Class used for template storage and atoms distance calculation
	"""

	def __init__(self,filename):
		atm = compile(r"^ATOM.{9}CA.{7}(?P<resid>.{4})(?P<x>.{12})(?P<y>.{8})(?P<z>.{8})")
		f = open(filename)
		for line in f.readlines():
			data = atm.match(line) # parse pdb file
			if data:
				print data.groupdict(),data.group('resid')
		f.close()




a = Template("./playground/2pcy.pdb")
