function [dout]=LO_bpsosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);

d_bp=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);

d_bpsosvmf = LO_sosvmf(d_bp,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);

dout=d_bpsosvmf-LO_fk_dip(d_bpsosvmf,w);
end