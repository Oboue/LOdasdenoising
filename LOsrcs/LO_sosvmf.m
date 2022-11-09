function [dout] = LO_sosvmf(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth)
% LO_sosvmf: structure oriented space-varying median filter (SOSVMF) method 
%
% Oboue et al., 2022
%
% INPUT
% din: input data (nt*nx)
% niter: number of nonlinear iterations
% liter: number of linear iterations (in divn)
% order: accuracy order
% eps_dv: eps for divn  (default: 0.01)
% eps_cg: eps for CG    (default: 1)
% tol_cg: tolerence for CG (default: 0.000001)
% rect:  smoothing radius (ndim*1)
% verb: verbosity flag

% adj:      adjoint flag
% add:      adding flag
% dip: slope (2D array)
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
%
%
% OUTPUT
% dout:     output data

dipn=str_dip2d(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb);
[~,ds]=LO_pwsmooth_lop_mf(adj,add,dipn,[],n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,din,[]);
dout=reshape(ds,n1,n2);

return
