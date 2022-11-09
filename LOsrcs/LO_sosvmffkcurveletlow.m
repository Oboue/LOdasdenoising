function [dout]=LO_sosvmffkcurveletlow(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,niter2,eps1,verb2);

d0=LO_sosvmffkcurvelet(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);

nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb2);

end