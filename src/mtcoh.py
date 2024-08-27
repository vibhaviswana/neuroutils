#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computes either multi-tapered phase locking value (PLV) or coherence between 
a given pair of signals. 
Copyright 2019-24 Vibha Viswanathan. All rights reserved.

INPUTS:
x: Signal 1 (size: number of trials x number of time points)
y: Signal 2 (size: number of trials x number of time points)
nw: Time-half-bandwidth product (scalar)
fs: Sampling frequency in Hz (scalar)
doPLV: Indicates whether to compute PLV (doPLV == true) or coherence (doPLV == false)
fmax: Maximum frequency of interest in Hz (scalar)

OUTPUTS:
outp: PLV or coherence (size: 1 x number of frequency points)
freqs: Frequencies (size: 1 x number of frequency points)

REFERENCES:
1. Hannan EJ (1970). Inference about spectra. In: Multiple time series,
Vol 1, pp 245–324. Hoboken, NJ: Wiley.
2. Thomson D (1982) Spectrum estimation and harmonic analysis. Proc IEEE 70:1055–1096.
3. Lachaux, J., Rodriguez, E., Martinerie, J., and Varela, F. (1999). “Measuring 
phase synchrony in brain signals,” Hum. Brain Mapp. 8(4), 194–208.
4. Zhu, L., Bharadwaj, H., Xia, J., and Shinn-Cunningham, B. (2013). “A comparison 
of spectral magnitude and phase-locking value analyses of the frequency-following
response to complex tones,” J. Acoust. Soc. Am. 134(1), 384–395.
5. Dobie RA, Wilson MJ (1994) Objective detection of 40 Hz auditory evoked 
potentials: phase coherence vs. magnitude-squared co-herence. Electroencephalogr
Clin Neurophysiol 92:405–413.
6. Slepian D (1978) Prolate spheroidal wave functions, Fourier analysis, and 
uncertainty V: the discrete case. Bell Syst Tech J 57:1371– 1430.
"""

import numpy as np
from anlffr.dpss import dpss_windows
from scipy.fft import fft

def mtcoh(x, y, nw, fs, doPLV, fmax):
    ntrials = x.shape[0]
    ntime = x.shape[1]
    ntapers = 2*nw-1
    [list_tapers,temp] = dpss_windows(ntime,nw,ntapers)
    nfft = int(2**np.ceil(np.log2(ntime)))
    freqs = np.arange(0,nfft)*fs/nfft
    freqinds = (freqs<=fmax)
    freqs = freqs[freqinds]
    nfreqs = freqs.shape[0]
    PLV = np.zeros((ntapers,nfreqs))
    Sxy = np.zeros((ntapers,nfreqs))
    Sxx = np.zeros((ntapers,nfreqs))
    Syy = np.zeros((ntapers,nfreqs))
    
    for k in np.arange(ntapers):
        tap = np.tile(list_tapers[k,:],(ntrials,1))
        Xf = fft(tap*x,n=nfft,axis=-1)
        Yf = fft(tap*y,n=nfft,axis=-1)
        Xf = Xf[:,freqinds]
        Yf = Yf[:,freqinds]
        if np.logical_not(doPLV):
            Sxy[k,:] = np.abs(np.mean(Xf*np.conjugate(Yf),axis=0))
            Sxx[k,:] = np.abs(np.mean(Xf*np.conjugate(Xf),axis=0))
            Syy[k,:] = np.abs(np.mean(Yf*np.conjugate(Yf),axis=0)) 
        elif doPLV:
            PLV[k,:] = np.abs(np.mean((Xf/np.abs(Xf))*(np.conjugate(Yf)/np.abs(Yf)), axis=0)) 
        #end  
    #end
    
    if np.logical_not(doPLV):
        outp = np.mean(Sxy,axis=0)/np.sqrt(np.mean(Sxx,axis=0)*np.mean(Syy,axis=0))
    elif doPLV:
        outp = np.mean(PLV,axis=0)
    #end

    return [outp,freqs]
#end

