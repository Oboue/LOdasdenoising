function [dout]=LO_bpsosvmffklow(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,rec,niter2,eps1,verb2)

d_bp=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);

d_bpsosvmf = LO_sosvmf(d_bp,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);

d0=d_bpsosvmf-LO_fk_dip(d_bpsosvmf,w);

%%
nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb2);

end