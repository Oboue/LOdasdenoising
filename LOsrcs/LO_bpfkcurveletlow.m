function [dout]=LO_bpfkcurveletlow(din,dt,flo,fhi,nplo,nphi,phase,verb0,w,n1,n2,c1,c2,c3,niter1,rec,niter2,eps1,verb2);
% LO_bpfkcurveletlow: BP+FK+Curvelet+LOW (LO) method

% Oboue et al., 2022
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
% w:        half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% d_est:    estimated signal using another denoising method. d_est corresponds to the input noisy data when the curvelet method is used as a single denosing step.
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% n1:       first dimension
% n2:       second dimension
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration

% rec:      3-D vector denoting smooth radius
% niter2:   number of CG iterations
% eps1:     regularization parameter, default 0.0
% verb2:    verbosity flag (default: 0)

% OUTPUT
% dout:     output data
%% BP
d_bp=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
%% BP+FK
d_fk=d_bp-LO_fk_dip(d_bp,w);
d_est=d_fk;
%% BP+FK+curvelet
d0=LO_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1);
%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb2);
end 