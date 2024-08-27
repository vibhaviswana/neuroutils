#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Test mtcoh.py 
# Copyright 2019-24 Vibha Viswanathan. All rights reserved.

import numpy as np
import pylab as pl
from mtcoh import mtcoh

fs = 1000.
ntrials = 100
t = np.arange(0., 1., 1./fs)
x = np.random.randn(ntrials, t.shape[0])

rhos = np.asarray([0, 0.2, 0.4, 0.8, 1])

pl.figure()

for rho in rhos:
    rhoflip  = (1 - rho ** 2.) ** 0.5
    y = rho * x + rhoflip * np.random.randn(ntrials, t.shape[0])
    
    nw = 3
    doPLV = False
    fmin, fmax = 0., 200.
    C, f = mtcoh(x, y, nw, fs, doPLV, fmin, fmax)
    
    pl.plot(f, C)
#end