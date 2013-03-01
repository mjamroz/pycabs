#!/usr/bin/env python
# 2013, Michal Jamroz, public domain. http://biocomp.chem.uw.edu.pl

import os, random, pylab, glob, pycabs, numpy as np, multiprocessing as mp
# first of all, download pyCABS and set self.FF = "" to the FF directory with cabs files
# to compile CABS, use: g77 -O2 -ffloat-store -static -o cabs CABS.f 
# to compile lattice model builder, use g77 -O2 -ffloat-store -static build_cabs61.f

sequence, secstr = pycabs.parseDSSPOutput("1bnr.dssp") # define file with secondary structure definition of define sequence and secondary structure in sequence,secstr variables respectively
name = "barnase" # name for the project. Script will create subdirectories with this name as suffix
template = ["../1bnr.pdb"] # set path to the start structure (here - native structure). If user want to start from random chain, set template=[]. Note that path is relative to simulation directory
independent_runs=5  # set number of independent simulations for each temperature
temp_from = 1.0     # define range of simulation temperatures, here is 1.0 - 2.8 with interval of 0.1
temp_to  = 2.8
temp_interval = 0.1
temperatures=np.arange(temp_from,temp_to,temp_interval)

def runCABS(temperature):
   global name, sequence,secstr,template,independent_runs
   here = os.getcwd()
   for i in range(independent_runs):
      temp = "%5.3f" %(temperature)
      dir_name= name+"_"+str(i)+"_T"+temp  # create unique name for simulation dir
      a = pycabs.CABS(sequence,secstr,template,dir_name) 
      a.rng_seed = random.randint(1,10000) # set random generator seed for each independent simulation
      a.createLatticeReplicas(replicas=1)  # create lattice model for CABS
      a.modeling(Ltemp=temperature,Htemp=temperature, phot=100,cycles=100) # start modeling. phot is CABS microcycle, cycles variable is CABS macrocycle (how often write to the trajectory file)
      os.chdir(here)

pool = mp.Pool()  # it use all available CPUs on workstation. If user want to use only - for example two - CPUs, set pool = mp.Pool(2)
pool.map(runCABS,temperatures) # run simulations in parallel way, each simulation on each available CPU

# postprocessing (comment out two lines above to avoid starting over simulations. If you want to only plot with other labels, etc. )

cv = np.empty([independent_runs,len(temperatures)])
avgene = np.empty([independent_runs,len(temperatures)])
for j in range(independent_runs):
   for i in range(len(temperatures)):
      t = temperatures[i]
      temp = "%5.3f" %(t)
      e_path = os.path.join(name+'_'+str(j)+'_T'+temp,'ENERGY') # path constructed in the same way like dir_name in runCABS definition above
      energy = np.fromfile(e_path,sep='\n') # read CABS energies for each trajectory model
      cv[j][i] = np.std(energy)  # calculate standard deviation (numpy std function) of energy, for each independent simulation
      avgene[j][i] = np.mean(energy) # calculate mean (numpy mean function) of energy


mean_sigma = np.mean(cv,axis=0)    # average over independent simulations
stddev_sigma = np.std(cv,axis=0)
mean_ene = np.mean(avgene,axis=0)  # calculate of standard deviations and mean values over independent simulations
stddev_ene = np.std(avgene,axis=0)

# plotting data with pylab python module. Read matplotlib manual (http://matplotlib.org/) for explanation of below code
pylab.ylabel(r'Standard deviation of energy' )
pylab.xlabel(r'Temperature, $T$')
pylab.xlim(temp_from,temp_to)
for i in range(independent_runs):
    pylab.plot(temperatures, cv[i], '.')
pylab.errorbar(temperatures,mean_sigma,yerr=stddev_sigma,fmt='o-')
pylab.savefig("stdE_barnase.png")
pylab.close()

pylab.ylabel(r'Mean energy' )
pylab.xlabel(r'Temperature, $T$')
pylab.xlim(temp_from,temp_to)
for i in range(independent_runs):
    pylab.plot(temperatures, avgene[i], '.')
pylab.errorbar(temperatures,mean_ene,yerr=stddev_ene,fmt='o-')
pylab.savefig("meanE_barnase.png")
pylab.close()



# as a results user will get a lot of barnase* subdirectories, stdE_barnase.png and meanE_barnase.png files.
# EOF
