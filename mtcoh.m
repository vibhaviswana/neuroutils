function [outp,freqs] = mtcoh(x, y, nw, fs, doPLV, fmax)
%
% Computes either multi-tapered phase locking value (PLV) or coherence between 
% a given pair of signals. 
% Copyright 2019-23 Vibha Viswanathan. All rights reserved.
%
% INPUTS:
% x: Signal 1 (size: number of trials x number of time points)
% y: Signal 2 (size: number of trials x number of time points)
% nw: Time-half-bandwidth product (scalar)
% fs: Sampling frequency in Hz (scalar)
% doPLV: Indicates whether to compute PLV (doPLV == true) or coherence (doPLV == false)
% fmax: Maximum frequency of interest in Hz (scalar)
%
% OUTPUTS:
% outp: PLV or coherence (size: 1 x number of frequency points)
% freqs: Frequencies (size: 1 x number of frequency points)
%
% REFERENCES:
% 1. Hannan EJ (1970). Inference about spectra. In: Multiple time series,
% Vol 1, pp 245–324. Hoboken, NJ: Wiley.
% 2. Thomson D (1982) Spectrum estimation and harmonic analysis. Proc IEEE 70:1055–1096.
% 3. Lachaux, J., Rodriguez, E., Martinerie, J., and Varela, F. (1999). “Measuring 
% phase synchrony in brain signals,” Hum. Brain Mapp. 8(4), 194–208.
% 4. Zhu, L., Bharadwaj, H., Xia, J., and Shinn-Cunningham, B. (2013). “A comparison 
% of spectral magnitude and phase-locking value analyses of the frequency-following
% response to complex tones,” J. Acoust. Soc. Am. 134(1), 384–395.
% 5. Dobie RA, Wilson MJ (1994) Objective detection of 40 Hz auditory evoked 
% potentials: phase coherence vs. magnitude-squared co-herence. Electroencephalogr
% Clin Neurophysiol 92:405–413.
% 6. Slepian D (1978) Prolate spheroidal wave functions, Fourier analysis, and 
% uncertainty V: the discrete case. Bell Syst Tech J 57:1371– 1430.

ntrials = size(x,1);
ntime = size(x,2);
ntaps = 2*nw-1;
list_taps = dpss(ntime,nw,ntaps);
nfft = 2^nextpow2(ntime);
freqs = (0:(nfft-1))*fs/nfft;
freqs = freqs(freqs<=fmax);
nfreqs = numel(freqs);
coh = zeros(ntaps,nfreqs);

for k = 1:ntaps
    tap = (repmat(list_taps(:,k),1,ntrials))';
    Xf = fft(tap.*x,nfft,2);
    Yf = fft(tap.*y,nfft,2);
    Xf = Xf(:,1:nfreqs);
    Yf = Yf(:,1:nfreqs);
    Yf = conj(Yf);
    if ~doPLV
        coh(k,:) = abs(mean(Xf.*Yf,1)./mean(abs(Xf).*abs(Yf),1)); 
    elseif doPLV
        coh(k,:) = abs(mean((Xf./abs(Xf)) .* (Yf./abs(Yf)),1)); 
    end
    
end
outp = mean(coh,1);


