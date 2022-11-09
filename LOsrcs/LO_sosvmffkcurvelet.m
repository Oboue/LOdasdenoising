function [dout]=LO_sosvmffkcurvelet(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1)
% LO_sosvmffkcurvelet: SOSVMF+FK+Curvelet method 
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

% w: half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

% d_est:    estimated signal using another denoising method. d_est corresponds to the input noisy data when the curvelet method is used as a single denosing step.
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% n1:       first dimension
% n2:       second dimension
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration
%% SOSVMF
d_sosvmf=LO_sosvmf(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
%% SOSVMF+FK
d_sosvmffk=d_sosvmf-LO_fk_dip(d_sosvmf,w);
%% SOSVMF+FK+Curvelet
d_sosvmffkcurvelet=LO_curvelet(din,d_sosvmffk,n1,n2,c1,c2,c3,niter1);
dout=d_sosvmffkcurvelet;
end