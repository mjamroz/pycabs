#!/usr/bin/env python
from pylab import *
import time

ion()

y = zeros(1)		# energy
x = zeros(1)
line, = plot(x,y,'--')
xlabel('CABS time step')
ylabel('CABS energy')

file = open('/tmp/energy', 'r')
i=0
while 1:
	where = file.tell()
	l = file.readline()
	if not l:
		file.seek(where)
		time.sleep(1)
	else:

		y = np.append(y,float(l))
		x = np.append(x,i)
		axis([0, amax(x)+1, amin(y)-5, amax(y)+5 ])

		i+=1
        line.set_ydata(y)  # update the data
        line.set_xdata(x)
        draw()
        
