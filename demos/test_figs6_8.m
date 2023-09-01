% Demo for an advanced median filter (AMF) for improving the signal-to-noise ratio of seismological datasets
% 

% Script to plot Figures 6, 7 and 8


%  Copyright (C) Oboue et al., 2023
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
%% BP+SOMF+FK
tic
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

LO_snr(d,dn) % Noisy
% LO_snr(d,d_dct) % 
% LO_snr(d,d_sgk) % 
LO_snr(d,d1_bpmffk) % 
LO_snr(d,d_bpsomffk) % 
load d_afk % load afk data 
LO_snr(d,d_afk) % 
LO_snr(d,d5) % 

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
text(n2/0.6,-0.07,'Clean','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Clean','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,3,2);LO_imagesc(dn,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'Noisy','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,3);LO_imagesc(d1_bpmffk,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,4);LO_imagesc(d_bpsomffk,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,5);LO_imagesc(d_afk,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'AFK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,3,6);LO_imagesc(d5,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

% dn1=dn-d_dct;
% dn2=dn-d_sgk;
dn3=dn-d1_bpmffk;
dn4=dn-d_bpsomffk;
dn5=dn-d_afk;
dn6=dn-d5;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,4,1);LO_imagesc(dn3,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,2);LO_imagesc(dn4,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,3);LO_imagesc(dn5,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'AFK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

subplot(2,4,4);LO_imagesc(dn6,95,1,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.75+ngap+n2,-0.04,'Raw','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2+ngap*2+n2*2,-0.04,'Removed noise','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(0.1,0.95,labels{5},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');

rect=[10,10,1];niter=20;eps=0;verb=0;
% [simid1]=LO_localsimi(dn-d_dct,d_dct,rect,niter,eps,verb);
% [simid2]=LO_localsimi(dn-d_sgk,d_sgk,rect,niter,eps,verb);
[simid3]=LO_localsimi(dn-d1_bpmffk,d1_bpmffk,rect,niter,eps,verb);
[simid4]=LO_localsimi(dn-d_bpsomffk,d_bpsomffk,rect,niter,eps,verb);
[simid5]=LO_localsimi(dn-d_afk,d_afk,rect,niter,eps,verb);
[simid6]=LO_localsimi(dn-d5,d5,rect,niter,eps,verb);

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
ax(1)=subplot(2,4,5); imagesc(x,t,simid3); hold on;
text(n2/0.60,-0.07,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(e)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colormap(ax(1),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

% subplot(2,4,6);LO_imagesc(simid4,95,1,x,t);colormap(jet);caxis([0 0.75]);
ax(2)=subplot(2,4,6); imagesc(x,t,simid4); hold on;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.60,-0.07,'BP+SOMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(f)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colormap(ax(2),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

ax(3)=subplot(2,4,7); imagesc(x,t,simid5); hold on;
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.60,-0.07,'AFK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(g)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colormap(ax(3),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');

ax(4)=subplot(2,4,8); imagesc(x,t,simid6); hold on;
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.60,-0.07,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(h)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colormap(ax(4),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');


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
%% AFK

trace3 = d_afk(:,inds);
%figure;
% subplot(4,1,3);
plot(t,trace3,'cyan','linewidth',4); hold on
% legend('BP+SOMF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',2,'Fontsize',40,'Fontweight','bold');

%% Proposed LO method
% inds=30;
trace4 = d5(:,inds);
% subplot(4,1,4)
plot(t,trace4,'r','linewidth',4); hold on
legend('Clean','Raw','BP+MF+FK','BP+SOMF+FK','AFK','LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-0.5 0.5])
set(gca,'Linewidth',6,'Fontsize',40,'Fontweight','bold');
annotation(gcf,'line',[0.17109375 0.145703125000001],...
    [0.867807153965785 0.867807153965785],...
    'Color',[1 0.411764705882353 0.16078431372549],...
    'LineWidth',10);
annotation(gcf,'textbox',...
    [0.169359375 0.846811819595645 0.16071875 0.0482115085536547],...
    'String','P wave arrival',...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'line',[0.30625 0.302734375],...
    [0.897356143079316 0.150077760497667],...
    'Color',[1 0.411764705882353 0.16078431372549],...
    'LineWidth',10);

% computation time size (3001*334)
% BP+MF+FK = 2.362394 seconds.
% BP+SOMF+FK =  350.766550 seconds.
% Proposed LO = 518.928361 seconds.
% 
%  
% size(d5)