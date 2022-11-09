function [dout]=LO_fk_dip_low(din,w,rec,niter2,eps1,verb);

[d0]=din-LO_fk_dip(din,w);

%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb);
end
