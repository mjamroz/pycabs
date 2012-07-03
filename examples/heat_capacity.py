#!/usr/bin/env python

import pycabs
import multiprocessing as mp
import os
def runCABS(temperature):
	global name, sequence,secstr,template
	
	# function for running CABS with different temperatures
	# it will compute in directory name+_+temperature
	here = os.getcwd()
	a = pycabs.CABS(sequence,secstr,template,name+"_"+str(temperature))
	a.createLatticeReplicas()
	a.modeling(Ltemp=temperature,Htemp=temperature, phot=1,cycles=1)
	os.chdir(here)


# lets create thread pool 

name = "fnord"
template=["/home/hydek/pycabs/playground/2pcy_CA.pdb"]
sss =  pycabs.parsePorterOutput("/home/hydek/pycabs/proba/playground/porter.ss") # suppose we have porter prediction of sec. str.
sequence = sss[0]
secstr = sss[1]

# creating working pool 
pool = mp.Pool(processes=2)
temperatures=xrange(1,5) # ranges of temperature

pool.map(runCABS,temperatures)

# wait ... and process name*/ENERGY files
