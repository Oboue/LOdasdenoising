function [dout]= LO_curveletlow(din,n1,n2,c1,c2,c3,niter1,rec,niter2,eps1,verb);
% LO_curveletlow: Curvelet + LOW (LO) method 
%
%Oboue et al., 2022
%
%INPUT

% din:      input data
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
% verb1:    verbosity flag (default: 0)

% OUTPUT
% dout:     output data

%% Curvelet
d_est=din;
d0=LO_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1);

%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb);
end