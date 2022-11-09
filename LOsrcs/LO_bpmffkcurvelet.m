function [dout]=LO_bpmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,nfw,ifb,axis,w,n1,n2,c1,c2,c3,niter1);

% LO_bpmffkcurvelet: BP+MF+FK+Curvelet method.

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
% nfw:      window size 
% ifb:      if use padded boundary (if not, zero will be padded)
% axis:     temporal sampling interval
% w:        half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% d_est:    estimated signal using another denoising method. d_est corresponds to the input noisy data when the curvelet method is used as a single denosing step.
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% n1:       first dimension
% n2:       second dimension
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration

% OUTPUT
% dout:     output data

%%
% BP
d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
% BP+MF
d2=LO_mf(d1,nfw,ifb,axis);
%BP+MF+FK
d3=d2-LO_fk_dip(d2,w);
%%BP+MF+FK+Curvelet
dout=LO_curvelet(din,d3,n1,n2,c1,c2,c3,niter1);
end