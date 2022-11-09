function [dout]=LO_sosvmflow(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,rec,eps1,niter2,verb2)
% LO_sosvmflow: SOSVMF + LOW (LO) method 
%
% Oboue et al., 2022
%
% INPUT
% din:      input data (nt*nx)
% niter:    number of nonlinear iterations
% liter:    number of linear iterations (in divn)
% order:    accuracy order
% eps_dv:   eps for divn  (default: 0.01)
% eps_cg:   eps for CG    (default: 1)
% tol_cg:   tolerence for CG (default: 0.000001)
% rect:     smoothing radius (ndim*1)
% verb:     verbosity flag

% adj:      adjoint flag
% add:      adding flag
% dip:      slope (2D array)
% w1:       weight
% n1:       trace length
% n2:       number of traces
% ns:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);
% 
% ndn: size of dn (n1*n2)
% nds: size of ds (n1*n2)
% 
% type_mf=0 (MF) or 1 (SVMF) 
% ifsmooth=1 (if smooth) or 0 (only MF)
% 
% dn: model  (1D array) noisy data
% ds: data   (1D array)smoothed data

% rec:      3-D vector denoting smooth radius
% niter2:   number of CG iterations
% eps1:     regularization parameter, default 0.0
% verb1:    verbosity flag (default: 0)

% OUTPUT
% dout:     output data

%% SOSVMF 
  d0 = LO_sosvmf(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);

%% Local orthogonalization operation 
  nois_0=din-d0; % compute the initial noise section
 [dout,nois2,low]=LO_localortho(d0,nois_0,rec,eps1,niter2,verb2);
end