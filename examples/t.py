#!/usr/bin/env python

import time,os,threading

class Monitor(threading.Thread):
	def __init__(self,filename):
		threading.Thread.__init__(self)
		self.daemon = False #: if True, it will terminate when script terminates		
		self.f = filename
		self.first_mtime = os.stat(filename).st_mtime
		self.kill = False
		fb = open(filename)
		self.buf = [i.strip() for i in fb.readlines()] # get whole file
		self.position = fb.tell() # and save position of eof
		fb.close()
		
	def terminate(self):
		self.kill=True
		
	def run(self):
		print self.buf
		self.buf=[]
		while 1:
			if self.kill:
				return
			time.sleep(1)
			mtime = os.stat(self.f).st_mtime
			if mtime!=self.first_mtime:
				self.first_mtime = mtime
				fb = open(self.f)
				fb.seek(self.position)
				#print self.position
				while 1:
					l = fb.readline().strip()
					if(l!=""):
						self.buf.append(l)
					else:
						self.position = fb.tell()
						break
				if "END" in self.buf: # reach end of file
					break
				else:
					print self.buf
					self.buf=[]
					
						
						
		
m=Monitor("/tmp/energy")
m.start()

#m.terminate()
