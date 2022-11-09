function [dout]=LO_curvelet(din,d_est,n1,n2,c1,c2,c3,niter1)
% LO_curvelet: curvelet method
%
%
%
%INPUT

% din:      input data
% d_est:    estimated signal using another denoising method. d_est corresponds to the input noisy data when the curvelet method is used as a single denosing step.
% c1:       Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
% n1:       first dimension
% n2:       second dimension
% c2:       Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
% c3:       Thresholding parameter (alpha)
% niter1:   Number of iteration

% OUTPUT
% dout:     output data

% Modify by Oboue et al., 2022

F=ones(n1,n2);                                   
X=fftshift(ifft2(F))*sqrt(prod(size(F)));  
C=LO_fdct_wrapping(X,0,c2);
% Compute norm of curvelets (exact)
E=cell(size(C));
for s=1:length(C)
    E{s}=cell(size(C{s}));
    for w=1:length(C{s})
         A=C{s}{w};
         E{s}{w}=sqrt(sum(sum(A.*conj(A)))/prod(size(A)));    
    end
end   
Cdn=LO_fdct_wrapping(din,c1,c2);     
Smax=length(Cdn);
Sigma0=c3*median(median(abs(Cdn{Smax}{1})))/0.58;    
Sigma=Sigma0;
sigma=[Sigma,linspace(2.5*Sigma,0.5*Sigma,niter1)];
Sigma=sigma(1);
Cdn=LO_fdct_wrapping(d_est,c1,c2);   
  
    Ct=Cdn;
    for s=2:length(Cdn)
        thresh=Sigma+Sigma*s;
        for w=1:length(Cdn{s})
            Ct{s}{w}=Cdn{s}{w}.*(abs(Cdn{s}{w})>thresh*E{s}{w});
        end
    end  
dout=real(LO_ifdct_wrapping(Ct,c1,n1,n2));
end