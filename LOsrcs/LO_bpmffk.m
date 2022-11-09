function [dout]=LO_bpmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,nfw,ifb,axis,w);

d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);

d_mf=LO_mf(d1,nfw,ifb,axis);
dout=d_mf-LO_fk_dip(d_mf,w);
end 