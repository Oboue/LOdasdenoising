function [dout]=LO_bpmffklow(din,dt,flo,fhi,nplo,nphi,phase,verb0,nfw,ifb,axi,w,rec,niter2,eps1,verb1);

d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);

d_mf=LO_mf(d1,nfw,ifb,axi);
d0=d_mf-LO_fk_dip(d_mf,w);
%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb1);
end