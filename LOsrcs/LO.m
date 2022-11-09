function [dout] = LO(din,par);
%  LO: Local orthogonalization method (LO)
%
% By Oboue et al., 2022
%
% INPUT
% din:      input data
% par:      parameter file
% IN par:
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
% w:        half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration
%
% rec:      3-D vector denoting smooth radius 
% eps1:     regularization parameter, default 0.0
% niter2:   number of CG iterations
% verb:     verbosity flag (default: 0) 
%
% OUTPUT
% dout:     output data
%%
[n1,n2]=size(din);

if isfield(par,'dt')
    dt=par.dt;
else
    dt=0.004;
end

if isfield(par,'flo')
    flo=par.flo;
else
    flo=0;
end

if isfield(par,'fhi')
    fhi=par.fhi;
else
    fhi=1/dt/2;
end

if isfield(par,'nplo')
    nplo=par.nplo;
else
    nplo=6;
end

if isfield(par,'nphi')
    nphi=par.nphi;
else
    nphi=6;
end

if isfield(par,'phase')
    phase=par.phase;
else
    phase=0;
end

if isfield(par,'verb0')
    verb0=par.verb0;
else
    verb0=0;
end

if isfield(par,'niter')
    niter=par.niter;
else
    niter=5;
end

if isfield(par,'liter')
    liter=par.liter;
else
    liter=20;
end

if isfield(par,'order1')
    order1=par.order1;
else
    order1=2;
end

if isfield(par,'eps_dv')
    eps_dv=par.eps_dv;
else
    eps_dv=0.01;
end

if isfield(par,'eps_cg')
    eps_cg=par.eps_cg;
else
    eps_cg=1;
end

if isfield(par,'tol_cg')
    tol_cg=par.tol_cg;
else
    tol_cg=0.000001;
end

if isfield(par,'rect')
    rect=par.rect;
else
    rect=25;
end

if isfield(par,'verb1')
    verb1=par.verb1;
else
    verb1=1;
end

if isfield(par,'adj')
    adj=par.adj;
else
    adj=0;
end

if isfield(par,'add')
    add=par.add;
else
    add=0;
end

if isfield(par,'n1')
    n1=par.n1;
end

if isfield(par,'n1')
    n2=par.n2;
end

if isfield(par,'ns')
    ns=par.ns;
else
    ns=8;
end

if isfield(par,'order2')
    order2=par.order2;
else
    order2=2;
end

if isfield(par,'eps')
    eps=par.eps;
else
    eps=0.01;
end

if isfield(par,'ndn')
    ndn=par.ndn;
else
    ndn=n1*n2;
end

if isfield(par,'nds')
    nds=par.nds;
else
    nds=n1*n2;
end

if isfield(par,'type_mf')
    type_mf=par.type_mf;
else
    type_mf=1;
end

if isfield(par,'ifsmooth')
    ifsmooth=par.ifsmooth;
else
    ifsmooth=0;
end

if isfield(par,'w')
    w=par.w;
else
    w=0;
end

if isfield(par,'c1')
    c1=par.c1;
else
    c1=1;
end

if isfield(par,'c2')
    c2=par.c2;
else
    c2=1;
end

if isfield(par,'c3')
    c3=par.c3;
else
    c3=0.5;
end

if isfield(par,'niter1')
    niter1=par.niter1;
else
    niter1=10;
end

if isfield(par,'rec')
    rec=par.rec;
else
    rec=10;
end

if isfield(par,'eps1')
    eps1=par.eps1;
else
    eps1=0;
end

if isfield(par,'niter2')
    niter2=par.niter2;
else
    niter2=1;
end

if isfield(par,'verb')
    verb=par.verb;
else
    verb=1;
end

[d0]=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);

% LO_bandpasssosvmffkcurvelet: robust denoising method composed of the bandpass filter, SOSVMF, FK, and the curvelet methods. 
% This function estimates the initially denoised signal d0

%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb);

% LO_localortho: LO_localortho compensates for the signal leakage energy
% causes by ''LO_bandpasssosvmffkcurvelet'' method.
% dout is the final result from the LO method.
return
