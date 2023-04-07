% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figure 3

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.

%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%-------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../LO_data/');
addpath('../LOsrcs/');
addpath('../seistr/');
%% Process the DAS seismic data with strong signal energy corrupted by strong noise using the LO method

NOs=[1,20,10,25,11,2];
labels={...                                             %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
    'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
    'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*3+2*ngap;

for ii=2
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
eq1=d1;
%% Denosing using the BP+MF+FK method 
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d1=LO_bandpass(eq,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
ns=8;  % spray radius
%
% tic
d_mf=LO_mf(d1,ns*2+1,1,2);
% toc

%FK
d1_bpmffk1=d_mf-LO_fk_dip(d_mf,0.02);%
toc
%%
% d_bp=dl_bandpass(eq1,0.0005,0,200,6,6,0,0);%
% % d_bpfk=d_bp-dl_fk_dip(d_bp,0.02);%
% d_bpfk=d_bp;
% d1=d_bpfk;
% figure;dl_imagesc([eq,d,eq-d]);

% load(strcat('/Users/chenyk/dasdenoising/mat_bpsomffk/eq-',num2str(ii),'.mat'));
% load data/real11.mat
% d=d1;
% figure;dl_imagesc([eq,d,eq-d]);

% dn=dn(:,100:800);
% d=d(:,100:800);
%% denoising
%% patch size l1*l2
% l1=32;l2=1;
% c1=l1;c2=64;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% %% DCT dictionary (dctmtx will generates orthogonal transform)
% dct=zeros(c1,c2);
% for k=0:1:c2-1
%     V=cos([0:1:c1-1]'*k*pi/c2);
%     if k>0
%         V=V-mean(V);
%     end
%     dct(:,k+1)=V/norm(V);
% end
% D=dct;
% % DCT=kron(dct,dct);%2D DCT dictionary (64,256)
% %% plot the first 64 atoms
% % figure;
% % for ia=1:16
% %     subplot(4,4,ia);plot(dct(:,ia));
% % end
% 
% % decompose the image into patches:
% X=dl_patch(d1,1,l1,1,l1/2,1);
% 
% % OMP using DCT
% nd=size(X,2);
% K=3;
% ph=2;
% % tic
% for i2=1:nd
%     G(:,i2)=dl_omp0(D,X(:,i2),K);
% end
% % toc
% 
% %further constrain it to be sparser
% G=dl_pthresh(G,'ph',ph);
% X2=D*G;
% 
% [n1,n2]=size(d1);
% d_dct1=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
% %% SGK
% param.T=K;      %sparsity level
% param.D=D;      %initial D
% param.niter=30; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% param.K=c2;     %number of atoms, dictionary size
% % tic
% [Dsgk,Gsgk]=dl_sgk(X,param); 
% % toc
% Gsgk0=Gsgk;
% Gsgk=dl_pthresh(Gsgk0,'ph',ph);
% X11=Dsgk*Gsgk;
% d_sgk1=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
% tic
d1=LO_bandpass(eq,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
% ns=8;  % spray radius
% %
% tic
% d_mf=LO_mf(d1,ns*2+1,1,2);
% toc

% BP+SOMF 

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=0;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
% tic
d1=LO_somf(d1,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% toc

% BP+SOMF+FK
% tic
d_bpsomffk1=d1-LO_fk_dip(d1,0.02);%
toc
% dn11=eq1-d_bpsomffk_2;
%% Denosing using the BP method
%  Parameter tuning for the BP method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
% tic
% d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
%%  Denosing using the BP+SOSVMF method
% Parameter tuning for the BP method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
% tic
% d2=LO_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.04;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% tic
% d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
% toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.05;             % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
% tic
% d4=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% toc
%
%
par.dt=dt;                  % sampling
par.flo=flo;                % Low frequency in band, default is 0
par.fhi=fhi;                % High frequency in band, default is Nyquist
par.nplo=nplo;              % number of poles for low cutoff
par.nphi=nphi;              % number of poles for high cutoff
par.phase=phase;            % y: minimum phase, n: zero phase
par.verb0=verb0;            % verbosity flag
%
par.niter=niter;            % number of nonlinear iterations
par.liter=liter;            % liter: number of linear iterations (in divn)
par.order1=order1;          % order: accuracy order
par.eps_dv=eps_dv;          % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;          % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;          % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=rect(1);        % rect:  smoothing radius (ndim*1)
par.rect(2)=rect(2);            % "      "        "
par.rect(3)=rect(3);        % "      "        "
par.verb1=1;            % verbosity flag

par.adj=adj;                % adjoint flag
par.add=add;                % adding flag
par.ns=ns;                  % spray radius
par.order2=order2;          % PWD order
par.eps=eps;                % regularization (default:0.01);
par.ndn=n1*n2;              % size of dn (n1*n2)
par.nds=n1*n2;              % size of ds (n1*n2)
par.type_mf=type_mf;        % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;      % 1 (if smooth) or 0 (only MF);
% 
par.w=0.05;                    % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                  % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c1;                  % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.08;                  % Thresholding parameter (alpha)
par.niter1=niter1;           % Number of iteration
%
rec = zeros(3, 1);          % 3-D vector denoting smooth radius 
par.rec(1) = 500;
par.rec(2) = 500;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb=1;                 % verbosity flag (default: 0) 
%
tic
d_LO1=LO(eq1,par);
% d5=d_LO1;
toc
% 
end
%% Raw data
t=[0:n1-1]*0.0005;
inds=20;
trace = eq1(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); hold on
plot(t,trace,'black','linewidth',3); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk1(:,inds);
%figure;
plot(t,trace1,'b','linewidth',3); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk1(:,inds);
%figure;
plot(t,trace2,'g','linewidth',3); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO1(:,inds);
plot(t,trace3,'r','linewidth',3); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','LO');
% xlim([1 2500])
xlabel('Times (s)')
ylabel('Amplitude')
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(1,1,1);LO_imagesc(d_LO1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.03,'Raw','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
%
% a = 13.1857; b = 14.1726; c = 15.18

% comp1=[eq1,zeros(n1,ngap),d_LO1,zeros(n1,ngap),eq1-d_LO1]; 
%% plot the waveforms weaker events

t=[0:n1]*0.0005;
ngap=50;
x=1:n2*3+2*ngap;
figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(2,3,1);LO_imagesc(eq1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.03,'Raw','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
%% Local similarity
rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(eq1-d_dct1,d_dct1,rect,niter,eps,verb);
% [simi2]=LO_localsimi(eq1-d_sgk1,d_sgk1,rect,niter,eps,verb);
[simi3]=LO_localsimi(eq1-d1_bpmffk1,d1_bpmffk1,rect,niter,eps,verb);
[simi4]=LO_localsimi(eq1-d_bpsomffk1,d_bpsomffk1,rect,niter,eps,verb);
[simi5]=LO_localsimi(eq1-d_LO1,d_LO1,rect,niter,eps,verb);
% 
%% Process the DAS seismic data with weak events 
%
eq=zeros(2000,960);
[n1,n2]=size(eq);
 
for ii=5;
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
eq2=d1;

%% Denosing using the BP+MF+FK method 
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d1=LO_bandpass(eq2,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
ns=8;  % spray radius
%
tic
d_mf=LO_mf(d1,ns*2+1,1,2);
toc

%FK
d1_bpmffk2=d_mf-LO_fk_dip(d_mf,0.02);%
%%
% d_bp=dl_bandpass(eq2,0.0005,0,200,6,6,0,0);%
% % d_bpfk=d_bp-dl_fk_dip(d_bp,0.02);%
% d_bpfk=d_bp;
% d1=d_bpfk;
% % figure;dl_imagesc([eq,d,eq-d]);
% %% denoising
% %% patch size l1*l2
% l1=32;l2=1;
% c1=l1;c2=64;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% %% DCT dictionary (dctmtx will generates orthogonal transform)
% dct=zeros(c1,c2);
% for k=0:1:c2-1
%     V=cos([0:1:c1-1]'*k*pi/c2);
%     if k>0
%         V=V-mean(V);
%     end
%     dct(:,k+1)=V/norm(V);
% end
% D=dct;
% DCT=kron(dct,dct);%2D DCT dictionary (64,256)
%% plot the first 64 atoms
% figure;
% for ia=1:16
%     subplot(4,4,ia);plot(dct(:,ia));
% end

% decompose the image into patches:
% X=dl_patch(d1,1,l1,1,l1/2,1);
% 
% % OMP using DCT
% nd=size(X,2);
% K=3;
% ph=2;
% tic
% for i2=1:nd
%     G(:,i2)=dl_omp0(D,X(:,i2),K);
% end
% toc
% 
% %further constrain it to be sparser
% G=dl_pthresh(G,'ph',ph);
% X2=D*G;
% 
% [n1,n2]=size(d1);
% d_dct2=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
% %% SGK
% param.T=K;      %sparsity level
% param.D=D;      %initial D
% param.niter=30; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% param.K=c2;     %number of atoms, dictionary size
% tic
% [Dsgk,Gsgk]=dl_sgk(X,param); 
% toc
% Gsgk0=Gsgk;
% Gsgk=dl_pthresh(Gsgk0,'ph',ph);
% X11=Dsgk*Gsgk;
% d_sgk2=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);

%% BP+SOMF+FK
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d2=LO_bandpass(eq2,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%

% BP+SOMF 

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=0;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
tic
d2=LO_somf(d2,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc

% BP+SOMF+FK
tic
d_bpsomffk2=d2-LO_fk_dip(d2,0.02);%
toc
% dn22=eq2-d_bpsomffk2;
%%
%% Denosing using the BP method
%  Parameter tuning for the BP method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
% tic
% d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
%%  Denosing using the BP+SOSVMF method
% Parameter tuning for the BP method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=15;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
% tic
% d2=LO_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.04;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% tic
% d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
% toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.05;             % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
% tic
% d4=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% toc
%
%
par.dt=dt;                  % sampling
par.flo=flo;                % Low frequency in band, default is 0
par.fhi=fhi;                % High frequency in band, default is Nyquist
par.nplo=nplo;              % number of poles for low cutoff
par.nphi=nphi;              % number of poles for high cutoff
par.phase=phase;            % y: minimum phase, n: zero phase
par.verb0=verb0;            % verbosity flag
%
par.niter=niter;            % number of nonlinear iterations
par.liter=liter;            % liter: number of linear iterations (in divn)
par.order1=order1;          % order: accuracy order
par.eps_dv=eps_dv;          % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;          % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;          % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=rect(1);        % rect:  smoothing radius (ndim*1)
rect(2)=rect(2);            % "      "        "
par.rect(3)=rect(3);        % "      "        "
par.verb1=verb1;            % verbosity flag

par.adj=adj;                % adjoint flag
par.add=add;                % adding flag
par.ns=10;                  % spray radius
par.order2=order2;          % PWD order
par.eps=eps;                % regularization (default:0.01);
par.ndn=n1*n2;              % size of dn (n1*n2)
par.nds=n1*n2;              % size of ds (n1*n2)
par.type_mf=type_mf;        % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;      % 1 (if smooth) or 0 (only MF);
% 
par.w=0.08;                    % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                  % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c1;                  % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                  % Thresholding parameter (alpha)
par.niter1=niter1;          % Number of iteration

rec = zeros(3,1);          % 3-D vector denoting smooth radius 
par.rec(1) = 500;
par.rec(2) = 500;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb =verb;                % verbosity flag (default: 0) 

%
tic
d_LO2=LO(eq2,par);
d5=d_LO2;
toc
% 
end

t=[0:n1-1]*0.0005;
inds=20;
trace = eq2(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); hold on
plot(t,trace,'black','linewidth',4); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk2(:,inds);
%figure;
plot(t,trace1,'b','linewidth',4); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk2(:,inds);
%figure;
plot(t,trace2,'g','linewidth',4); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO2(:,inds);
plot(t,trace3,'r','linewidth',4); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','BP+SOSVMF+FK+Curvelet+LOW (LO)');
% xlim([1 2500])
xlabel('Times (s)')
ylabel('Amplitude')
set(gca,'Linewidth',2,'Fontsize',35,'Fontweight','bold');

%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))

% 13.6925 14.8476 15.72
%% plot the waveforms weaker events

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% subplot(2,3,1);LO_imagesc(eq2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'Raw','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

% subplot(2,3,2);LO_imagesc(d_dct2,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.04,'DCT','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

% subplot(2,3,2);LO_imagesc(d_sgk2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,3);LO_imagesc(d1_bpmffk2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'BP+MF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,4);LO_imagesc(d_bpsomffk2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'BP+SOMF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,5);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

% dn12=eq2-d_dct2;
% dn22=eq2-d_sgk2;
dn32=eq2-d1_bpmffk2;
dn42=eq2-d_bpsomffk2;
dn52=eq2-d_LO2;

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% % subplot(2,5,1);LO_imagesc(dn12,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.65,-0.04,'DCT','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,1);LO_imagesc(dn22,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'SGK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,2);LO_imagesc(dn32,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,3);LO_imagesc(dn42,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,4);LO_imagesc(dn52,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi12]=LO_localsimi(eq2-d_dct2,d_dct2,rect,niter,eps,verb);
% [simi22]=LO_localsimi(eq2-d_sgk2,d_sgk2,rect,niter,eps,verb);
[simi32]=LO_localsimi(eq2-d1_bpmffk2,d1_bpmffk2,rect,niter,eps,verb);
[simi42]=LO_localsimi(eq2-d_bpsomffk2,d_bpsomffk2,rect,niter,eps,verb);
[simi52]=LO_localsimi(eq2-d_LO2,d_LO2,rect,niter,eps,verb);

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% % subplot(2,5,6);LO_imagesc(simi12,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.65,-0.02,'DCT','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
% % text(-275,-0.04,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,5);LO_imagesc(simi22,95,2,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'SGK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.05,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,6);LO_imagesc(simi32,95,2,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.05,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,7);LO_imagesc(simi42,95,2,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.05,'(h)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,8);LO_imagesc(simi52,95,2,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.05,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%% Process the DAS seismic data with intermediate signal energy corrupted by strong noise using the LO method
eq=zeros(2000,960);
[n1,n2]=size(eq);

for ii=6;
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
eq3=d1;

%% Denosing using the BP+MF+FK method 
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag

%
tic
d1=LO_bandpass(eq3,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%
ns=8;  % spray radius
%
tic
d_mf=LO_mf(d1,ns*2+1,1,2);
toc

%FK
d1_bpmffk3=d_mf-LO_fk_dip(d_mf,0.02);%
%%
% d_bp=dl_bandpass(eq3,0.0005,0,200,6,6,0,0);%
% % d_bpfk=d_bp-dl_fk_dip(d_bp,0.02);%
% d_bpfk=d_bp;
% d1=d_bpfk;
% % figure;dl_imagesc([eq,d,eq-d]);
% 
% % load(strcat('/Users/chenyk/dasdenoising/mat_bpsomffk/eq-',num2str(ii),'.mat'));
% % load data/real11.mat
% % d=d1;
% % figure;dl_imagesc([eq,d,eq-d]);
% % dn=dn(:,100:800);
% % d=d(:,100:800);
% %% denoising
% %% patch size l1*l2
% l1=32;l2=1;
% c1=l1;c2=64;%size of the 1D cosine dictionary (if c2>c1, overcomplete)
% %% DCT dictionary (dctmtx will generates orthogonal transform)
% dct=zeros(c1,c2);
% for k=0:1:c2-1
%     V=cos([0:1:c1-1]'*k*pi/c2);
%     if k>0
%         V=V-mean(V);
%     end
%     dct(:,k+1)=V/norm(V);
% end
% D=dct;
% % DCT=kron(dct,dct);%2D DCT dictionary (64,256)
% %% plot the first 64 atoms
% % figure;
% % for ia=1:16
% %     subplot(4,4,ia);plot(dct(:,ia));
% % end
% 
% % decompose the image into patches:
% X=dl_patch(d1,1,l1,1,l1/2,1);
% 
% % OMP using DCT
% nd=size(X,2);
% K=3;
% ph=2;
% tic
% for i2=1:nd
%     G(:,i2)=dl_omp0(D,X(:,i2),K);
% end
% toc
% 
% %further constrain it to be sparser
% G=dl_pthresh(G,'ph',ph);
% X2=D*G;
% 
% [n1,n2]=size(d1);
% d_dct3=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
% %% SGK
% param.T=K;      %sparsity level
% param.D=D;      %initial D
% param.niter=30; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% param.K=c2;     %number of atoms, dictionary size
% tic
% [Dsgk,Gsgk]=dl_sgk(X,param); 
% toc
% Gsgk0=Gsgk;
% Gsgk=dl_pthresh(Gsgk0,'ph',ph);
% X11=Dsgk*Gsgk;
% d_sgk3=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);
%% BP
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d3=LO_bandpass(eq3,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%

% ns=8;  % spray radius
% %
% tic
% d_mf=LO_mf(d1,ns*2+1,1,2);
% toc

% BP+SOMF 

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

ns=8;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=0;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
tic
d3=LO_somf(d3,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc

% BP+SOMF+FK
tic
d_bpsomffk3=d3-LO_fk_dip(d3,0.02);%
toc

% dn33=eq3-d_bpsomffk_3;

%% Denosing using the LO method 
%
%% Denosing using the BP method
%  Parameter tuning for the BP method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
% tic
% d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
%%  Denosing using the BP+SOSVMF method
% Parameter tuning for the BP method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=50;                   % rect:  smoothing radius (ndim*1)
rect(2)=50;                   % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=15;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
% tic
% d2=LO_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.04;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% tic
% d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
% toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.05;             % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

% 
% tic
% d4=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% toc
%
%%
par.dt=dt;                  % sampling
par.flo=flo;                % Low frequency in band, default is 0
par.fhi=fhi;                % High frequency in band, default is Nyquist
par.nplo=nplo;              % number of poles for low cutoff
par.nphi=nphi;              % number of poles for high cutoff
par.phase=phase;            % y: minimum phase, n: zero phase
par.verb0=verb0;            % verbosity flag
%
par.niter=niter;            % number of nonlinear iterations
par.liter=liter;            % liter: number of linear iterations (in divn)
par.order1=order1;          % order: accuracy order
par.eps_dv=eps_dv;          % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=eps_cg;          % eps_cg: eps for CG    (default: 1)
par.tol_cg=tol_cg;          % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=rect(1);        % rect:  smoothing radius (ndim*1)
rect(2)=rect(2);            % "      "        "
par.rect(3)=rect(3);        % "      "        "
par.verb1=verb1;            % verbosity flag

par.adj=adj;                % adjoint flag
par.add=add;                % adding flag
par.ns=8;                  % spray radius
par.order2=order2;          % PWD order
par.eps=eps;                % regularization (default:0.01);
par.ndn=n1*n2;              % size of dn (n1*n2)
par.nds=n1*n2;              % size of ds (n1*n2)
par.type_mf=type_mf;        % 0 (MF) or 1 (SVMF)
par.ifsmooth=ifsmooth;      % 1 (if smooth) or 0 (only MF);
% 
par.w=w;                    % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                  % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c1;                  % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                  % Thresholding parameter (alpha)
par.niter1=niter1;           % Number of iteration
%
rec = zeros(3,1);          % 3-D vector denoting smooth radius 
par.rec(1) = 250;
par.rec(2) = 250;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb =1;                 % verbosity flag (default: 0) 
%
tic
d_LO3=LO(eq3,par);
% d5=d_LO3
toc
% 
end

t=[0:n1-1]*0.0005;
inds=20;
trace = eq3(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); hold on
plot(t,trace,'black','linewidth',4); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk3(:,inds);
%figure;
plot(t,trace1,'b','linewidth',4); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk3(:,inds);
%figure;
plot(t,trace2,'g','linewidth',4); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO3(:,inds);
plot(t,trace3,'r','linewidth',4); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','BP+SOSVMF+FK+Curvelet+LOW (LO)');
% xlim([1 2500])
xlabel('Times (s)')
ylabel('Amplitude')
set(gca,'Linewidth',2,'Fontsize',35,'Fontweight','bold');

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(1,1,1);LO_imagesc(d_LO3,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.03,'Raw','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

%% S/N estimation
a=20*log(var(trace1))/log(var(trace));
b=20*log(var(trace2))/log(var(trace))
c=20*log(var(trace3))/log(var(trace))

% a = 12.0894, b = 12.8922 c = 15.1902
%%
% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
%     subplot(2,3,1);LO_imagesc(eq3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     ylabel('Time (s)','Fontsize',12,'fontweight','bold');
%     xlabel('Channel','Fontsize',12,'fontweight','bold');
%     set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     text(n2/0.65,-0.03,'Raw','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     text(-200,-0.05,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
%     
%     % subplot(2,3,2);LO_imagesc(d_dct3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
%     % xlabel('Channel','Fontsize',16,'fontweight','bold');
%     % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     % text(n2/0.65,-0.04,'DCT','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     % text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     
%     subplot(2,3,2);LO_imagesc(d_sgk3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     ylabel('Time (s)','Fontsize',12,'fontweight','bold');
%     xlabel('Channel','Fontsize',12,'fontweight','bold');
%     set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     text(n2/0.65,-0.03,'SGK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     text(-200,-0.05,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
%     
%     subplot(2,3,3);LO_imagesc(d1_bpmffk3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     ylabel('Time (s)','Fontsize',12,'fontweight','bold');
%     xlabel('Channel','Fontsize',12,'fontweight','bold');
%     set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     text(n2/0.65,-0.03,'BP+MF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     text(-200,-0.05,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
%     
%     subplot(2,3,4);LO_imagesc(d_bpsomffk3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     ylabel('Time (s)','Fontsize',12,'fontweight','bold');
%     xlabel('Channel','Fontsize',12,'fontweight','bold');
%     set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     text(n2/0.65,-0.03,'BP+SOMF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     text(-200,-0.05,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
%     
%     subplot(2,3,5);LO_imagesc(d_LO3,95,2,x,t);caxis([-25,25]); hold on
%     % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
%     ylabel('Time (s)','Fontsize',12,'fontweight','bold');
%     xlabel('Channel','Fontsize',12,'fontweight','bold');
%     set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
%     text(n2/0.65,-0.03,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%     text(-200,-0.05,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

% dn13=eq3-d_dct3;
% dn23=eq3-d_sgk3;
dn33=eq3-d1_bpmffk3;
dn43=eq3-d_bpsomffk3;
dn53=eq3-d_LO3;

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% % subplot(2,5,1);LO_imagesc(dn13,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.65,-0.02,'DCT','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.05,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,1);LO_imagesc(dn23,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'SGK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,2);LO_imagesc(dn33,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,3);LO_imagesc(dn43,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,4);LO_imagesc(dn53,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi13]=LO_localsimi(eq3-d_dct3,d_dct3,rect,niter,eps,verb);
% [simi23]=LO_localsimi(eq3-d_sgk3,d_sgk3,rect,niter,eps,verb);
[simi33]=LO_localsimi(eq3-d1_bpmffk3,d1_bpmffk3,rect,niter,eps,verb);
[simi43]=LO_localsimi(eq3-d_bpsomffk3,d_bpsomffk3,rect,niter,eps,verb);
[simi53]=LO_localsimi(eq3-d_LO3,d_LO3,rect,niter,eps,verb);

figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,5,6);LO_imagesc(simi13,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'DCT','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.04,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

% subplot(2,4,5);LO_imagesc(simi23,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.02,'SGK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.05,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,4,6);LO_imagesc(simi33,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.02,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-275,-0.05,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,4,7);LO_imagesc(simi43,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.02,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-275,-0.05,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,4,8);LO_imagesc(simi53,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.02,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-275,-0.05,'(h)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');


% dn3=eq3-d_LO3;
% % comp3=[eq3,zeros(n1,ngap),d_LO3,zeros(n1,ngap),eq3-d_LO3]; 
% %%
% %% Plot figures 
% 
% % figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% % subplot(3,1,1);LO_imagesc(comp2,95,1,x,t);
% % ylabel('Time (s)','Fontsize',13,'fontweight','bold');
% % xlabel('Channel','Fontsize',13,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
% % text(n2/2,-0.04,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap+n2,-0.04,'LO','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
% % annotation(gcf,'textarrow',[0.429894179894177 0.437830687830685],...
% %     [0.84830097087378 0.884708737864072],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % annotation(gcf,'textarrow',[0.479790662553815 0.476821912553815],...
% %     [0.800315260175598 0.772075761775026],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'String',{'  Visible','signal'},...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % subplot(3,1,2);LO_imagesc(comp3,95,1,x,t);
% % ylabel('Time (s)','Fontsize',13,'fontweight','bold');
% % xlabel('Channel','Fontsize',13,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% % text(n2/2,-0.04,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap+n2,-0.04,'LO','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{6},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
% % annotation(gcf,'textarrow',[0.541005291005289 0.518518518518515],...
% %     [0.576456310679611 0.605582524271838],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % annotation(gcf,'textarrow',[0.482436165199317 0.479467415199317],...
% %     [0.472645357262958 0.444405858862385],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'String',{'  Visible','signal'},...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % subplot(3,1,3);LO_imagesc(comp1,95,1,x,t);
% % ylabel('Time (s)','Fontsize',13,'fontweight','bold');
% % xlabel('Channel','Fontsize',13,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% % text(n2/2,-0.04,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap+n2,-0.04,'LO','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{2},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
% % annotation(gcf,'textarrow',[0.417989417989416 0.443121693121689],...
% %     [0.317961165048544 0.288834951456304],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % annotation(gcf,'textarrow',[0.447089947089945 0.472222222222218],...
% %     [0.172330097087378 0.143203883495139],'Color',[1 0 0],'TextColor',[1 0 0],...
% %     'String',{'Visible','signal'},...
% %     'LineWidth',2,...
% %     'HorizontalAlignment','center',...
% %     'FontWeight','bold',...
% %     'FontSize',14,...
% %     'FontName','Helvetica Neue');
% % 
% % print(gcf,'-depsc','-r300','fig3.eps');
% 
% %%   plot the waveforms weaker events
% 
% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% subplot(3,3,1);LO_imagesc(eq2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,4);LO_imagesc(d_bpsomffk_2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,5);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,4);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% %% plot the removed noise sections
% 
% subplot(3,3,7);LO_imagesc(dn22,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,8);LO_imagesc(dn2,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,7);LO_imagesc(dn1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% %% plot the local similarity metrics
% 
% rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(eq2-d2,d2,rect,niter,eps,verb);
% [simi2]=LO_localsimi(eq2-d_LO2,d_LO2,rect,niter,eps,verb);
% % [simi3]=LO_localsimi(eq1-d_LO1,d_LO1,rect,niter,eps,verb);
% 
% % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(3,4,8);LO_imagesc(simi1,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
% text(n2/0.7,-0.015,'LO','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.03,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,4,9);LO_imagesc(simi2,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
% text(n2/0.7,-0.015,'LO','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-275,-0.03,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,10);LO_imagesc(dn2,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% %%   plot the waveforms moderate events
% 
% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% subplot(3,3,1);LO_imagesc(eq3,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,4);LO_imagesc(d_bpsomffk_3,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,5);LO_imagesc(d_LO3,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,4);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% %% plot the removed noise sections
% 
% subplot(3,3,7);LO_imagesc(dn33,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,8);LO_imagesc(dn3,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,7);LO_imagesc(dn1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% %%   plot the waveforms stronger events
% 
% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
% subplot(3,3,1);LO_imagesc(eq1,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,4);LO_imagesc(d_bpsomffk_1,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,5);LO_imagesc(d_LO1,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,4);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% %% plot the removed noise sections
% 
% subplot(3,3,7);LO_imagesc(dn11,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,8);LO_imagesc(dn1,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% xlabel('Channel','Fontsize',16,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(3,4,7);LO_imagesc(dn1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% 
% 
% 
% % % figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% % subplot(3,4,5);LO_imagesc(eq3,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % subplot(3,3,5);LO_imagesc(d_LO3,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % subplot(3,3,6);LO_imagesc(dn3,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % %%
% % subplot(3,3,7);LO_imagesc(eq1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % subplot(3,3,8);LO_imagesc(d_LO1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(h)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % subplot(3,3,9);LO_imagesc(dn1,95,2,x,t);caxis([-25,25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% % text(n2/0.7,-0.04,'Removed noise','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.1,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % %% Local similarity metric
% % 
% % rect=[30,30,1];niter=20;eps=0;verb=0;
% % [simi1]=LO_localsimi(eq2-d_LO2,d_LO2,rect,niter,eps,verb);
% % [simi2]=LO_localsimi(eq3-d_LO3,d_LO3,rect,niter,eps,verb);
% % [simi3]=LO_localsimi(eq1-d_LO1,d_LO1,rect,niter,eps,verb);
% % 
% % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% % subplot(1,3,1);LO_imagesc(simi1,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
% % text(n2/0.7,-0.015,'LO','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(-275,-0.03,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % subplot(1,3,2);LO_imagesc(simi2,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
% % text(n2/0.7,-0.015,'LO','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(-275,-0.03,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 
% % 
% % subplot(1,3,3);LO_imagesc(simi3,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',16,'fontweight','bold');
% % xlabel('Channel','Fontsize',16,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',16,'Fontweight','bold');
% % text(n2/0.7,-0.015,'LO','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(-275,-0.03,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % 

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(3,4,1);LO_imagesc(eq1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(a)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');

subplot(3,4,2);LO_imagesc(d1_bpmffk1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(b)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.364444444444444 0.366666666666667],...
    [0.772413793103445 0.739901477832502],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,3);LO_imagesc(d_bpsomffk1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(c)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.568888888888889 0.571111111111111],...
    [0.771428571428571 0.738916256157628],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,4);LO_imagesc(d_LO1,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(d)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{2},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.784444444444444 0.78],...
    [0.783251231527091 0.741871921182261],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,5);LO_imagesc(eq3,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{6},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');


subplot(3,4,6);LO_imagesc(d1_bpmffk3,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(f)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{6},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.407777777777778 0.407777777777778],...
    [0.557635467980295 0.533004926108365],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.368888888888889 0.388888888888889],...
    [0.458128078817734 0.44334975369457],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,7);LO_imagesc(d_bpsomffk3,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(g)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{6},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.571111111111111 0.591111111111111],...
    [0.454187192118225 0.439408866995061],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.618888888888889 0.618888888888889],...
    [0.558620689655171 0.533990147783241],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,8);LO_imagesc(d_LO3,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(h)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{6},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.824444444444444 0.824444444444444],...
    [0.560591133004923 0.535960591132993],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.778888888888889 0.798888888888889],...
    [0.456157635467978 0.441379310344814],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(3,4,9);LO_imagesc(eq2,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(i)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{5},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');

subplot(3,4,10);LO_imagesc(d1_bpmffk2,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(j)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{5},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.397777777777778 0.376666666666667],...
    [0.272906403940886 0.291625615763534],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.405555555555556 0.405555555555556],...
    [0.247290640394089 0.223645320197032],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.374444444444444 0.358888888888889],...
    [0.193103448275862 0.168472906403929],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,11);LO_imagesc(d_bpsomffk2,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(k)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{5},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.604444444444444 0.583333333333333],...
    [0.272906403940886 0.291625615763533],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.613333333333333 0.613333333333333],...
    [0.248275862068965 0.224630541871909],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.582222222222222 0.566666666666666],...
    [0.193103448275861 0.168472906403929],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);

subplot(3,4,12);LO_imagesc(d_LO2,95,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(l)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.95,labels{5},'color','b','Fontsize',5.75,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'arrow',[0.812222222222222 0.791111111111111],...
    [0.272906403940886 0.291625615763533],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.818888888888889 0.818888888888889],...
    [0.250246305418718 0.226600985221661],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
annotation(gcf,'arrow',[0.788888888888889 0.773333333333333],...
    [0.200985221674876 0.176354679802943],'Color',[1 0 0],'LineWidth',2,...
    'HeadWidth',25,...
    'HeadLength',15);
%%
% subplot(3,3,7);LO_imagesc(dn32,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.05,'(g)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,8);LO_imagesc(dn42,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.05,'(h)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(3,3,9);LO_imagesc(dn52,95,2,x,t);caxis([-25,25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.04,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.05,'(i)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
% %%
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(3,3,1);LO_imagesc(simi3,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+MF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(a)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,2);LO_imagesc(simi4,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+SOMF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(b)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,3);LO_imagesc(simi5,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'LO','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(c)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

% subplot(3,3,4);LO_imagesc(simi23,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.05,'(e)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,4);LO_imagesc(simi33,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+MF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(d)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,5);LO_imagesc(simi43,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+SOMF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(e)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,6);LO_imagesc(simi53,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'LO','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(f)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

% subplot(3,3,7);LO_imagesc(simi22,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.05,'(i)','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,7);LO_imagesc(simi32,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+MF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(g)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,8);LO_imagesc(simi42,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'BP+SOMF+FK','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(h)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

subplot(3,3,9);LO_imagesc(simi52,95,2,x,t);colormap(jet);caxis([0,0.5]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.06,'LO','color','k','Fontsize',22,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(i)','color','k','Fontsize',30,'fontweight','bold','HorizontalAlignment','center');

%%  Strong events
inds=20;
trace = eq1(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); hold on
subplot(3,1,1);
plot(t,trace,'black','linewidth',2); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk1(:,inds);
%figure;
plot(t,trace1,'b','linewidth',2); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk1(:,inds);
%figure;
plot(t,trace2,'g','linewidth',2); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO1(:,inds);
plot(t,trace3,'r','linewidth',2); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-100 100])
text(-0.08,125,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% BP+SOMF+FK method
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
%% Intermediate events
t=[0:n1-1]*0.0005;
inds=20;
trace = eq3(:,inds);
% length(trace)
subplot(3,1,2);
plot(t,trace,'black','linewidth',2); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk3(:,inds);
%figure;
plot(t,trace1,'b','linewidth',2); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk3(:,inds);
%figure;
plot(t,trace2,'g','linewidth',2); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO3(:,inds);
plot(t,trace3,'r','linewidth',2); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-100 100])
text(-0.08,125,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
%% Weak events

t=[0:n1-1]*0.0005;
inds=20;
trace = eq2(:,inds);
% length(trace)
subplot(3,1,3);
plot(t,trace,'black','linewidth',2); hold on
%% BP+SOMF+FK method
inds=20;
trace1 = d1_bpmffk2(:,inds);
%figure;
plot(t,trace1,'b','linewidth',2); hold on
%% BP+SOMF+FK method 
inds=20;
trace2 = d_bpsomffk2(:,inds);
%figure;
plot(t,trace2,'g','linewidth',2); hold on
%% Proposed LO method
inds=20;
trace3 = d_LO2(:,inds);
plot(t,trace3,'r','linewidth',2); hold on
legend('Raw','BP+MF+FK','BP+SOMF+FK','LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
set(gca,'Linewidth',2,'Fontsize',25,'Fontweight','bold');
ylim([-100 100])
 text(-0.08,125,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))

