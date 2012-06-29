#!/usr/bin/env python

import numpy as np

a = np.zeros(1)

for i in range(10):
	oldSize = a.shape[0]
	a.resize(oldSize+1)
	a[oldSize] = i
	print a
