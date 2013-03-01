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
	a.modeling(Ltemp=temperature,Htemp=temperature, phot=10,cycles=100)
	#remember to come back to `here` directory
	os.chdir(here)


#init these variables _before_ running cabs
name = "fnord"
# we have some template, it has to be as list
template=["/home/user/pycabs/playground/2pcy.pdb"] 
# suppose we have porter prediction of sec. str.
sss =  pycabs.parsePorterOutput("/home/user/pycabs/proba/playground/porter.ss") 
sequence = sss[0]
secstr = sss[1]
# now we have all data required to run CABS

temp_from = 1.5
temp_to  = 3.0
temp_interval = 0.1
temperatures=np.arange(temp_from,temp_to,temp_interval) # ranges of temperature

# create thread pool  with two parallel threads
pool = mp.Pool(processes=2)
pool.map(runCABS,temperatures) # run cabs threads

# HERE IS THE END OF PART WHERE WE RUN CABS in parallel fashion. 

# Now you can do something with output data, we'll calculate heat capacity, Cv:
cv = np.empty(len(temperatures))
for i in range(len(temperatures)):
	t = temperatures[i]
	e_path = os.path.join(name+'_'+str(t),'ENERGY')
	energy = np.fromfile(e_path,sep='\n') # read ENERGY data into array `energy`
	avg_energy2 = np.average(energy*energy) # <E^2>
	avg_energy = np.average(energy)		    # <E>^2
	cv[i] = np.std(energy)*np.std(energy)/(t*t) # (<E^2> - <E>^2) / T^2
# now we have heat capacity in cv array	

# ... and display plot
from pylab import *
xlabel(r'temperature $T$')
ylabel(r'heat capacity $C_v = (\left<E^2\right> - \left<E\right>^2)/T^2$' )
xlim(temp_from,temp_to) # xrange
plot(temperatures,cv)
show()

#remember that you have name+_+temperature directories, delete it or sth
