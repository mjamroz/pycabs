#!/usr/bin/env python
from pylab import *
from sys import argv
import time
import numpy as np
import pycabs

class Energy(pycabs.Calculate):
    def calculate(self,data):
        for i in data:
            self.out.append(float(i)) # ENERGY file contains one value in a row
            
out = []						
calc = Energy(out) # out is dynamically updated 
m=pycabs.Monitor(os.path.join(argv[1],"ENERGY"),calc)
m.daemon = True
m.start()



ion()
y = zeros(1)
x = zeros(1)
line, = plot(x,y)
xlabel('CABS time step')
ylabel('CABS energy')

while 1:
    time.sleep(1)
    y = np.asarray(out)
    x = xrange(0,len(out))
    axis([0, amax(x)+1, amin(y)-5, amax(y)+5 ])
    line.set_ydata(y)  # update the data
    line.set_xdata(x)
    draw()
