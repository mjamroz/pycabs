#!/usr/bin/env python
import pycabs
from Pycluster import *
from numpy import array,zeros

working_dir = "de_novo" # name of project 
templates = [] # deNOVO
a = pycabs.CABS('AAAAAAAAAA','HHHHHHHHHH',templates,working_dir) # initialize CABS, create required files
# DENOVO a.generateConstraints() 
a.createLatticeReplicas() # create start models from templates
a.modeling(Htemp=3.0,cycles=2,phot=1) # start modeling with default INP values and create TRAF.pdb when done
tr = a.getTraCoordinates() # load TRAF into memory and calculate RMSD all-vs-all : 


#calculating RMSD 2D array for clustering
distances = zeros((len(tr),len(tr)))
for i in range(len(tr)):
	for j in range(i,len(tr)):
		rms = pycabs.rmsd(tr[i],tr[j])
		distances[i][j] = distances[j][i] = rms
		
#save RMSD array as heat map
pycabs.heat_map(distances,"Protein model","Protein model","RMSD")		


# clustering by K-medoids method (with 5 clusters)
clusterid,error,nfound = kmedoids(distances,nclusters=5,npass=15,initialid=None)
print clusterid,error
clusterid,error,nfound = kcluster(distances,nclusters=5,npass=15)

# save cluster medoids to file
pycabs.saveMedoids(clusterid,a)
print clusterid,error
