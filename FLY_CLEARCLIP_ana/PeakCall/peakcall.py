#!/usr/bin/python
#this python script is for peak calling of CLASH data
#using cubic spline analysis
# usage: python peakcalling.py cluster.txt result.txt
#

import sys
import os
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter

f = open(sys.argv[1],"r")
local_max = []
for line in f:
	line = line.strip();
	cluster = line.split('\t')
	id = cluster[0]
	cst = cluster[1].split('|')
	coord = map(int,cst[0].split(','))
	height = map(int,cst[1].split(','))
	if max(height)<3:
		continue
	if len(coord)<5:
		continue
	yVar = gaussian_filter(height,1)
	cs = UnivariateSpline(coord, height,w=1.0/np.sqrt(yVar))
	x_base = np.arange(min(coord), max(coord)+1)
	y_interpolated = cs(x_base)
	y_prime = cs(x_base, 1)
	for i in range(len(y_prime)-1):
		left = y_prime[i]
		right = y_prime[i+1]
		if left > 0 and right < 0:
			cheight = y_interpolated[i]
			for j in range(len(coord)-1):
				if x_base[i]>= coord[j] and x_base[i]<=coord[j+1]:
					cheight = max(height[j-1],height[j],height[j+1])
			local_max.append((id,x_base[i], cheight,y_interpolated[i]))

w = open(sys.argv[2],'w')
for record in local_max:
	print >>w, record[0],"\t",record[1],"\t",record[2],"\t",record[3]
