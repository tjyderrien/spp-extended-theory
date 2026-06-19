#!/usr/bin/env python
#-*- coding: utf-8 -*-
from SimpleSPProutines import *
import pylab as plb

#reset

# 1. Build the 1D distribution of angular orientation

# Interface with ImageJ > OrientationJ ?  

#AngleData, AngleSpectrum

# 2. Get the tichkness of the distribution


# Read in data -- first 2 rows are header in this example. 

x = AngleData
y = AngleSpectrum

# Just an initial guess
mean = sum(x*y)
sigma = sum(y*(x - mean)**2)

# Fitting
def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
popt, pcov = curve_fit(gauss_function, x, y, p0 = [1, mean, sigma])

# 3. Adding a line into the experimental file. 