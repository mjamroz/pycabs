#!/usr/bin/env python
# 2013, Michal Jamroz, public domain. http://biocomp.chem.uw.edu.pl
import pycabs,os,numpy as np

name = "barnase" # set project name same like in folding_pathway.py script
max_sd_temperature=2.9 # read from the plot temperature for which maximum deviation of energy is observed
independent_runs=5     # set the same value like in folding_pathway.py script

# read trajectory of sidegroups (TRASG) from all independent simulations to the trajectory variable
trajectory = []
for j in range(independent_runs):
	temp = "%06.3f" %(max_sd_temperature)
	e_path = os.path.join(name+'_'+str(j)+'_T'+temp,'TRASG')
	trajectory += pycabs.loadSGCoordinates(e_path)[1000:] # second half of sidechains trajectory

contact = pycabs.contact_map(trajectory,7.0) # calculate averaged contact map over trajectories. Cutoff set to 7.0A.

# plot contact map with pylab. Read matplotlib manual (http://matplotlib.org/) for explanation of the code below
from pylab import xlabel,ylabel,pcolor,colorbar,savefig,xlim,ylim,cm
from numpy import indices
l=len(trajectory[0])
rows, cols = indices((l,l))
xlabel("Residue index")
xlim(0, len(contact))
ylim(0, len(contact))
ylabel("Residue index")
pcolor(contact, cmap=cm.gnuplot2_r,vmax=0.6) # vmax is range of colorbar. Here is set to 0.8, which gave white color for all values greater than 0.8
cb = colorbar()
cb.set_label("Fraction of contacts")
savefig("heatmap"+str(max_sd_temperature)+".png",dpi=600)

# optionally, write contact map values into the text file formatted for GNUplot
fw = open("contact_map.dat","w")
for i in range(len(contact)):
   for j in range(len(contact)):
      fw.write("%5d %5d %7.5f\n" %(i+1,j+1,contact[i][j]))
   fw.write("\n")
fw.close()

# example GNUplot script for plotting contact map of contact_map.dat file:
'''
set terminal unknown
plot 'contact_map.dat' using 1:2:3
set xrange[GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
set yrange[GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]


set terminal postscript eps enhanced color "Helvetica" 14
set output 'contact_map.eps'
set size ratio 1
unset key
set xlabel 'Residue index'
set ylabel 'Residue index'
set cbrange[:0.8]
set palette negative
plot 'contact_map.dat' with image
'''

# write it to the file.gp and run: gnuplot file.gp to get postscript file with heat map plot
