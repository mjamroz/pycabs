#!/usr/bin/env python

import time,os,threading

class Monitor(threading.Thread):
	"""
		Class for monitoring of CABS output data. You can run it and dynamically update output arrays with calculated results.
	"""	
	def __init__(self,filename,calculate):
		threading.Thread.__init__(self)
		self.daemon = False #: if True, it will terminate when script terminates		
		self.f = filename
		self.first_mtime = os.stat(filename).st_mtime
		self.kill = False
		self.fb = open(filename)
		self.buf = [i.strip() for i in self.fb] # get whole file
		self.position = self.fb.tell() # and save position of eof
		self.calc = calculate
		
		self.calc.calculate(self.buf) #update calculate object
		
	def terminate(self):
		self.kill=True
		
	def run(self):
		self.buf=[]
		while 1:
			if self.kill:
				return
			time.sleep(1)
			mtime = os.stat(self.f).st_mtime
			if mtime!=self.first_mtime:
				self.first_mtime = mtime
				self.fb.seek(self.position)
				for l in self.fb:
					self.buf.append(l.strip())
				self.position = self.fb.tell() #update EOF
				
				self.calc.calculate(self.buf) #update calculate object after each fortran FLUSH 
				self.buf=[]
					
class Calculate: # dziedzicz po tym rozne obliczenia
	def __init__(self,output):
		self.out = output
	def calculate(self,data):
		self.out.append('_'.join(data))	
								
		
		
		
