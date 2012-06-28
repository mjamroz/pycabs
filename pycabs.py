#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
		self.sequence = sequence
		self.ss = secondary_structure
		self.templates_fn = templates_filenames

	def p(self):
		print self.sequence


a = CABS("sdfsdf","sdfsdfsdf","sdf")
a.p()
