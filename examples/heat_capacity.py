#!/usr/bin/env python

import pycabs
import multiprocessing as mp

def runCABS(temperature,name,templates,sequence,ss):
	# function for running CABS with different temperatures
	# it will compute in directory name+_+temperature
	a = CABS(sequence,ss,templates,name+"_"+str(temperature))
	a.createLatticeReplicas()
	a.modeling(Ltemp=temperature,Htemp=temperature, phot=1,cycles=1)


# lets create thread pool 

name = "fnord"
temperatures=range(1,10)
pool = mp.Pool(processes=2)


# wait ...
