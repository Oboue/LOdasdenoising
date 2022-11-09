function [dout]=amf_bpfkcurveletlow(din,dt,flo,fhi,nplo,nphi,phase,verb0,w,n1,n2,c1,c2,c3,niter1,rec,niter2,eps1,verb);

d_bp=amf_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);

d_fk=amf_fk_dip(d_bp,w);
d_est=d_fk;

d0=amf_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1);
%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=amf_localortho(d0,nois_0,rec,niter2,eps1,verb);
end 