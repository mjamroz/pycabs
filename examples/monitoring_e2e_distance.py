#!/usr/bin/env python
from pylab import *
from sys import argv
import time
import numpy as np
import pycabs

class E2E(pycabs.Calculate):
    def calculate(self,data):
        models = self.processTrajectory(data)
        for m in models:
            first = m[0:3]
            last = m[-3:]

            x = first[0]-last[0]
            y = first[1]-last[1]
            z = first[2]-last[2]
            self.out.append(x*x+y*y+z*z)            
            
out = []						
calc = E2E(out) # out is dynamically updated 
m=pycabs.Monitor(argv[1]+"/TRAF",calc)
m.daemon = True
m.start()


ion()
y = zeros(1)
x = zeros(1)
line, = plot(x,y)
xlabel('CABS time step')
ylabel('square of end to end distance')

while 1:
    time.sleep(1)
    y = np.asarray(out)
    x = xrange(0,len(out))
    axis([0, amax(x)+1, amin(y)-5, amax(y)+5 ])
    line.set_ydata(y)  # update the data
    line.set_xdata(x)
    draw()
