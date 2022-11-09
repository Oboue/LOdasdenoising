function [dout]=LO_bandpasslow(din,dt,flo,fhi,nplo,nphi,phase,verb,rec,niter2,eps1,verb1);
% LO_bandpasslow: BP+LOW (LO) method

% Oboue et al., 2022

% Aug, 05, 2020
% by Yangkang Chen
%
% INPUT
% din:      input data
% dt:       sampling
% flo:      Low frequency in band, default is 0
% fhi:      High frequency in band, default is Nyquist
% nplo=6:   number of poles for low cutoff
% nphi=6:   number of poles for high cutoff
% phase=0:  y: minimum phase, n: zero phase
% verb=0:   verbosity flag

% rec:      3-D vector denoting smooth radius
% niter2:   number of CG iterations
% eps1:     regularization parameter, default 0.0
% verb1:    verbosity flag (default: 0)

% OUTPUT
% dout:     output data
%% BP
  d0=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb);
  %% Local orthogonalization operation 
  nois_0=din-d0; % compute the initial noise section
 [dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb1);
return
