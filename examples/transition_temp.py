#!/usr/bin/env python
import multiprocessing as mp
import os
import numpy as np

import pycabs

def runCABS(temperature):
	# global for simplify arguments 
	global name, sequence,secstr,template 
	
	# function for running CABS with different temperatures
	# it will compute in directory name+_+temperature
	here = os.getcwd() # since pycabs changing directories...
	a = pycabs.CABS(sequence,secstr,template,name+"_"+str(temperature))
	a.createLatticeReplicas(replicas=1)
	a.modeling(Ltemp=temperature,Htemp=temperature, phot=15,cycles=6)
	#remember to come back to `here` directory
	t = a.getTraCoordinates()
	os.chdir(here)
	return t

# define function for contacts 
# function return 1D array of 1/0
def contacts(model,cutoff=8.0, seq_sep=5):
	from numpy import sqrt
	contacts_arr=bytearray()
	prot_len = len(model)/3 # pycabs recognize model as 1D list of len. 3N
	for i in range(prot_len-seq_sep):
		for j in range(i+seq_sep,prot_len):
			x = model[3*i] - model[3*j]
			y = model[3*i+1] - model[3*j+1]
			z = model[3*i+2] - model[3*j+2]
			d = sqrt(x*x + y*y + z*z)
			b = 1 if d<cutoff else 0
			contacts_arr.append(b)
	return contacts_arr
def compare_contacts(contacts1,contacts2):
	s = 0.0
	for i in range(len(contacts1)):
		if contacts1[i]&contacts2[i]: # 0&0=0&1=0, 1&1=1
			s+=1
	return s	

#init these variables _before_ running cabs
name = "fnord"
# we have some template, it has to be as list
template=["/home/hydek/pycabs/playground/2pcy.pdb"] 
# suppose we have porter prediction of sec. str.
sss =  pycabs.parsePorterOutput("/home/hydek/pycabs/proba/playground/porter.ss") 
sequence = sss[0]
secstr = sss[1]
# now we have all data required to run CABS

temp_from = 0.1
temp_to  = 4.0
temp_interval = 0.15
temperatures=np.arange(temp_from,temp_to,temp_interval) # ranges of temperature

# create thread pool  with two parallel threads

pool = mp.Pool(processes=2)
trajectories = pool.map(runCABS,temperatures) # run cabs threads
	
# calculate contacts for native structure = "template" model for cabs
native = pycabs.parsePDBfile(template[0])
native_contacts = contacts(native)
nat_contacts_len = compare_contacts(native_contacts,native_contacts) # how many contacts in native structure
avg_nat_contacts = []
for trajectory in trajectories:
	nat_contacts_sum =0.0
	for model in trajectory:
		model_contacts = contacts(model)
		nat_contacts_sum += compare_contacts(native_contacts,model_contacts)
	avg_nat_contacts.append(nat_contacts_sum/(len(trajectory)*nat_contacts_len))



# display plot
from pylab import *
xlabel(r'temperature $T$')
ylabel(r'average number of native contacts, $<Q>$' )
xlim(temp_from,temp_to) # xrange
plot(temperatures,avg_nat_contacts)
show()
