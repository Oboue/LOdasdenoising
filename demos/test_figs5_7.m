% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets
% 
%  Copyright (C) Oboue et al., 2022
%-
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/

%  References

%  Oboue et al., 2022
%  Huang, G., M. Bai, Q. Zhao, W. Chen, and Y. Chen, 2021, Erratic noise suppression using iterative structure-oriented space-varying median filtering with sparsity constraint, Geophysical Prospecting, 69, 101-121.
%  Chen, Y., S. Zu, Y. Wang, and X. Chen, 2020, Deblending of simultaneous-source data using a structure-oriented space varying median filter, Geophysical Journal International, 222, 1805�1�723.
%  Zhao, Q., Q. Du, X. Gong, and Y. Chen, 2018, Signal-preserving erratic noise attenuation via iterative robust sparsity-promoting filter, IEEE Transactions on Geoscience and Remote Sensing, 56, 1558-0644.
%%
clc;clear;close all;
%% addpath
addpath('../LO_data/');
addpath('../LOsrcs/');
addpath('../seistr/');
%% load data clean data
% 
load micro_sf_3001_334_3.mat
d=LO_scale(data(:,:,1),2);
[n1,n2]=size(d);

% load horizontal + vertical + random + erratic noise
% load d_noise_micro.mat 
% load d_noise_micro.mat

% dt=0.0005;
% t=[0:nt-1]*dt;
% h=[0:nx-1]*dx+x0;
% ngap=50;
% x=1:n2*3+2*ngap;
% 
% load('fig2noise.mat')
% 
% randn('state',211111); % Horizontal
% nnn=LO_seisdither(nn,round(LO_meanf(20*randn(1,nx),20,1,2)));
% 
% dn=d+0.03*nnn;
% 
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,4,1);LO_imagesc(dn,95,1,x,t);
% 
% 
% randn('state',2111112); % Vertical 
% nn2=LO_seisdither(nn2,round(LO_meanfs(10*randn(1,nt),20,1,2,100)));
% 
% dn=dn + 0.2*nn2';
% 
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,4,1);LO_imagesc(dn,95,1,x,t);
% 
% 
% % Create the mask operator
% mask=rand(1,n2);
% mask(logical(mask<0.9))=0;
% mask(logical(mask>=0.9))=1;
% 
% % Adding erratic noise 
% 
% err_n=zeros(size(dn));
% for i=1:n1
%     randn('state',123456+i);
%     err_n(i,:)=0.2*randn(1,n2).*mask;
% end
% 
% % Adding random noise 
% 
% randn('state',201920);
% ran_n=0.07*randn(n1,n2);
% 
% dn=dn+err_n+ran_n;
% % din=dn;
% % 
% % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% % subplot(2,4,1);LO_imagesc(dn,95,1,x,t);
% % 
% LO_snr(d,dn) % Noisy
% % 
% save synthdas3.mat dn

load synthdas3.mat 
[n1,n2]=size(dn);
%%
% d_bp=dl_bandpass(dn,0.0005,0,200,6,6,0,0);%
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
% d_dct=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
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
% d_sgk=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);

%% Denosing using the BP+MF+FK method 
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
tic
d1=LO_bandpass(dn,dt,flo,fhi,nplo,nphi,phase,verb0);
ns=3;  % spray radius
d_mf=LO_mf(d1,ns*2+1,1,2);
d1_bpmffk=d_mf-LO_fk_dip(d_mf,0.05);%
toc
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
d1=LO_bandpass(dn,dt,flo,fhi,nplo,nphi,phase,verb0);
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
rect(1)=27;                   % rect:  smoothing radius (ndim*1)
rect(2)=7;                   % "      "        "
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
d_bpsomffk=d1-LO_fk_dip(d1,0.05);%
toc
%%  Denosing using the BP method
% Parameter tuning for the BP method
% tic
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
% tic
% d1=LO_bandpass(dn,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
%% Denosing using the BP+SOSVMF method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=20;                   % rect:  smoothing radius (ndim*1)
rect(2)=7;                    % "      "        "
rect(3)=1;                    % "      "        "
verb1=1;                      % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=5 ;                        % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)

% tic
% d2=LO_bandpasssosvmf(dn,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
% toc

%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.05;       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
% 
% tic
% d3=LO_bandpasssosvmffk(dn,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
% toc

% LO_snr(d,d2) % Noisy
% LO_snr(d,d3) % Noisy
% % 
%
% %% Denoising data using the BP+SOSVMF+FK+curvelet method 
% %  Parameter tuning：add the key parameters of the curvelet method
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=1.5;               % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
% tic
% d4=LO_bandpasssosvmffkcurvelet(dn,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% toc
%
%% Denosing using the LO method 
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
par.verb1=verb1;            % verbosity flag

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
par.w=w;                    % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=c1;                  % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=c1;                  % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=c3;                  % Thresholding parameter (alpha)
par.niter1=niter1;          % Number of iteration
%
rec = zeros(3, 1);          % 3-D vector denoting smooth radius 
par.rec(1) = 11;
par.rec(2) = 11;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb2=1;                % verbosity flag (default: 0) 
%
tic
d5=LO(dn,par);
toc
%%
load micro_sf_3001_334_3.mat
d=LO_scale(data(:,:,1),2);
[n1,n2]=size(d);
size(d)
% LO_snr(d,dn) % Noisy
% % LO_snr(d,d_dct) % 
% % LO_snr(d,d_sgk) % 
% LO_snr(d,d1_bpmffk) % 
% LO_snr(d,d_bpsomffk) % 
% LO_snr(d,d5) % 

dt=0.0005;
t=[0:nt-1]*dt;
h=[0:nx-1]*dx+x0;
ngap=50;
x=1:n2*3+2*ngap;

% comp1=[d,zeros(n1,ngap),dn]; 

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,3,1);LO_imagesc(d,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.6,-0.07,'Clean','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Clean','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,3,2);LO_imagesc(dn,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'Noisy','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,3);LO_imagesc(d1_bpmffk,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,4);LO_imagesc(d_bpsomffk,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,5);LO_imagesc(d5,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

% dn1=dn-d_dct;
% dn2=dn-d_sgk;
dn3=dn-d1_bpmffk;
dn4=dn-d_bpsomffk;
dn5=dn-d5;

% subplot(3,5,6);LO_imagesc(dn1,95,1,x,t);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Channel','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.60,-0.06,'DCT','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
% 
% subplot(3,5,7);LO_imagesc(dn2,95,1,x,t);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Channel','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.60,-0.06,'SGK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,3,1);LO_imagesc(dn3,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,2);LO_imagesc(dn4,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,3);LO_imagesc(dn5,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

rect=[10,10,1];niter=20;eps=0;verb=0;
% [simid1]=LO_localsimi(dn-d_dct,d_dct,rect,niter,eps,verb);
% [simid2]=LO_localsimi(dn-d_sgk,d_sgk,rect,niter,eps,verb);
[simid3]=LO_localsimi(dn-d1_bpmffk,d1_bpmffk,rect,niter,eps,verb);
[simid4]=LO_localsimi(dn-d_bpsomffk,d_bpsomffk,rect,niter,eps,verb);
[simid5]=LO_localsimi(dn-d5,d5,rect,niter,eps,verb);

%  figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(3,5,11);LO_imagesc(simid1,95,1,x,t);colormap(jet);caxis([0 1]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Channel','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',11,'Fontweight','bold');
% text(n2/0.60,-0.06,'DCT','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(k)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(3,5,12);LO_imagesc(simid2,95,1,x,t);colormap(jet);caxis([0 1]);
% ylabel('Time (s)','Fontsize',20,'fontweight','bold');
% xlabel('Channel','Fontsize',20,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',11,'Fontweight','bold');
% text(n2/0.60,-0.06,'SGK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% % text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(l)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,3,4);LO_imagesc(simid3,95,1,x,t);colormap(jet);caxis([0 0.75]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,5);LO_imagesc(simid4,95,1,x,t);colormap(jet);caxis([0 0.75]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOSVMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,6);LO_imagesc(simid5,95,1,x,t);colormap(jet);caxis([0 0.75]);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');


t=[0:n1-1]*0.0005;
inds=10;
trace = d(:,inds);
% length(d(:))
% length(trace)
figure('units','normalized','Position',[0.0 0.0 1 1],'color','w'); hold on
% subplot(4,1,1);
plot(t,trace,'g','linewidth',4); hold on
% legend('Raw');
% xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');

trace = dn(:,inds);
% length(d(:))
% length(trace)
% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); hold on
% subplot(4,1,1);
plot(t,trace,'black','linewidth',4); hold on
legend('Raw');
% xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');
%% BP+SOMF+FK method
t=[0:n1-1]*0.0005;
% inds=20;
trace1 = d1_bpmffk(:,inds);
%figure;
% subplot(4,1,2);
plot(t,trace1,'b','linewidth',4); hold on
% legend('BP+MF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');
%% BP+SOMF+FK method 
% inds=20;
trace2 = d_bpsomffk(:,inds);
%figure;
% subplot(4,1,3);
plot(t,trace2,'y','linewidth',4); hold on
% legend('BP+SOMF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');
%% Proposed LO method
% inds=30;
trace3 = d5(:,inds);
% subplot(4,1,4)
plot(t,trace3,'r','linewidth',4); hold on
legend('Clean','Raw','BP+MF+FK','BP+SOMF+FK','LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');


% computation time size (3001*334)
% BP+MF+FK = 2.362394 seconds.
% BP+SOMF+FK =  350.766550 seconds.
% Proposed LO = 518.928361 seconds.
% 
%  
size(d5)