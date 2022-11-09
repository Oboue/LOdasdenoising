function [dout] = LO_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% LO_bandpasssosvmf: Bandpass + SOSVMF method
%
% By Oboue et al., 2022
%
% INPUT
% din:      input data
% dt:       sampling
% flo:      Low frequency in band, default is 0
% fhi:      High frequency in band, default is Nyquist
% nplo:     number of poles for low cutoff
% nphi:     number of poles for high cutoff
% phase:    y: minimum phase, n: zero phase
% verb0:    verbosity flag
%
% niter:    number of nonlinear iterations
% liter:    number of linear iterations (in divn)
% order1:   accuracy order
% eps_dv:   eps for divn  (default: 0.01)
% eps_cg:   eps for CG    (default: 1)
% tol_cg:   tolerence for CG (default: 0.000001)
% rect:     smoothing radius (ndim*1)
% verb1:    verbosity flag
%
% adj:      adjoint flag
% add:      adding flag
% ns:       spray radius
% order2:   PWD order
% eps:      regularization (default:0.01);
% ndn:      size of dn (n1*n2)
% nds:      size of ds (n1*n2)
% type_mf:  0 (MF) or 1 (SVMF)
% ifsmooth: 1 (if smooth) or 0 (only MF)
%
%
% OUTPUT
% dout:     output data
%%
[d1]=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
% LO_bandpass: Bandpass filtering
%%
[dout] = LO_sosvmf(d1,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% LO_sosvmf: structure oriented space-varying median filter (SOSVMF) method 

return