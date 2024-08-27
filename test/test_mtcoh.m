
% Test mtcoh.m 
% Copyright 2019-24 Vibha Viswanathan. All rights reserved.

fs = 1000;
ntrials = 100;
t = (0:(1/fs):1);
x = randn(ntrials, numel(t));

rhos = [0, 0.2, 0.4, 0.8, 1];

figure()

for rho = rhos
    rhoflip  = (1 - rho^2)^0.5;
    y = rho * x + rhoflip * randn(ntrials, numel(t));
    
    nw = 3;
    doPLV = false;
    fmin = 0;
    fmax = 200;
    [C, f] = mtcoh(x, y, nw, fs, doPLV, fmin, fmax);
    
    plot(f, C);
    hold on;
end

