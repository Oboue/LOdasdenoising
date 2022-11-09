function [dout]=LO_fk_diplow(din,w,rec,niter2,eps1,verb);

% LO_fk_diplow: FK+ LOW (LO) method
%
%
% INPUT
% d: input data (2D)
% w: half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% rec:      3-D vector denoting smooth radius
% niter2:   number of CG iterations
% eps1:     regularization parameter, default 0.0
% verb1:    verbosity flag (default: 0)

% OUTPUT
% dout:     output data
%% FK 
[d0]=din-LO_fk_dip(din,w);

%% Local orthogonalization operation 
 nois_0=din-d0; % compute the initial noise section
[dout,nois2,low]=LO_localortho(d0,nois_0,rec,niter2,eps1,verb);
end
