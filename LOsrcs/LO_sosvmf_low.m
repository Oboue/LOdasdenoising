function [dout]=LO_sosvmf_low(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,rec,eps1,niter2,verb2)

  d0 = LO_sosvmf(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);

  %% Local orthogonalization operation 
  nois_0=din-d0; % compute the initial noise section
 [dout,nois2,low]=LO_localortho(d0,nois_0,rec,eps1,niter2,verb2);
end