% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figures 13, 14 and 15

%  Copyright (C) Oboue et al., 2023

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
%% Load the raw DAS seismic data corrupted by strong high-frequency noise, high-amplitude 

NOs=[1,20,10,25,11,2];
labels={...                                             %P-arrival sample NO from the SEGY file
    'FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy',... %24169
    'FORGE\_78-32\_iDASv3-P11\_UTC190426070723.sgy',... %24811
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062208.sgy',... %26090
    'FORGE\_78-32\_iDASv3-P11\_UTC190426110008.sgy',... %4921
    'FORGE\_78-32\_iDASv3-P11\_UTC190426062553.sgy',... %8934
    'FORGE\_78-32\_iDASv3-P11\_UTC190423182409.sgy'};   %4210
eq=zeros(2000,960);
[n1,n2]=size(eq)
t=[0:n1]*0.0005;
ngap0=1000;

ngap=50;
x=1:n2*3+2*ngap;
for ii=4
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
% din=d1;
% din=eq;

% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
% subplot(1,1,1);LO_imagesc(eq,100,2,x,t); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

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
% X=dl_patch(eq,1,l1,1,l1/2,1);
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
% [n1,n2]=size(eq);
% d_dct=dl_patch_inv(X2,1,n1,n2,l1,1,l1/2,1);
% 
% %% SGK
% param.T=K;      %sparsity level
% param.D=D;    %initial D
% param.niter=30; %number of K-SVD iterations to perform; default: 10
% param.mode=1;   %1: sparsity; 0: error
% param.K=c2;     %number of atoms, dictionary size
% % tic
% [Dsgk,Gsgk]=dl_sgk(X,param); 
% % toc
% Gsgk0=Gsgk;
% Gsgk=dl_pthresh(Gsgk0,'ph',ph);
% X11=Dsgk*Gsgk;
% d11=dl_patch_inv(X11,1,n1,n2,l1,1,l1/2,1);
% 
% inds=20:20:n2;
% traces=eq(:,inds);
% traces1=d11(:,inds);
% 
% dn=traces;
% d11trace=traces1;
% 
% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
% subplot(1,1,1);LO_imagesc(d11,100,2,x,t); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
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
d1=LO_bandpass(eq,dt,flo,fhi,nplo,nphi,phase,verb0);
% toc
%
ns=8;  % spray radius
%
% tic
d_mf=LO_mf(d1,ns*2+1,1,2);
% toc

%FK
d1_bpmffk=d_mf-LO_fk_dip(d_mf,0.02);%
toc

% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
% subplot(1,1,1);LO_imagesc(d1_bpmffk,100,2,x,t); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.65,-0.03,'SGK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

%% BP+SOMF+FK
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
d_bpsomffk=d1-LO_fk_dip(d1,0.02);%
toc

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(1,1,1);LO_imagesc(d_bpsomffk,100,2,x,t); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.03,'SGK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

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
%% 
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
% w=0.08;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
% tic
% d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
% toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.5;             % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
% tic
% d4=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
% toc
%
%% Denosing using the LO method 
%

par.dt=0.0005; % sampling
par.flo=0;     % Low frequency in band, default is 0
par.fhi=200;   % High frequency in band, default is Nyquist
par.nplo=6;    % number of poles for low cutoff
par.nphi=6;    % number of poles for high cutoff
par.phase=0;   % y: minimum phase, n: zero phase
par.verb0=0;   % verbosity flag
%
par.niter=2;                      % number of nonlinear iterations
par.liter=10;                     % liter: number of linear iterations (in divn)
par.order1=3;                     % order: accuracy order
par.eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
par.eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
par.tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
par.rect(1)=50;                   % rect:  smoothing radius (ndim*1)
par.rect(2)=50;                   % "      "        "
par.rect(3)=1;                    % "      "        "
par.verb1=1;                      % verbosity flag

par.adj=0;                        % adjoint flag
par.add=0;                        % adding flag
par.ns=12;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.03;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=.05;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 150;
par.rec(2) = 150;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;                 % verbosity flag (default: 0) 

%
tic
d_LO=LO(eq,par);
d4=d_LO;
toc

%% Raw data
t=[0:n1-1]*0.0005;
inds=90;
trace = eq(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); hold on
% subplot(5,1,1);
plot(t,trace,'black','linewidth',5); hold on
legend('Raw');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-90 90])
%  text(-0.08,125,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',30,'Fontweight','bold');
%% BP+SOMF+FK method
t=[0:n1-1]*0.0005;
% inds=20;
trace1 = d1_bpmffk(:,inds);
%figure;
% subplot(5,1,2);
plot(t,trace1,'b','linewidth',5); hold on
legend('BP+MF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-90 90])
% text(-0.08,600,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',30,'Fontweight','bold');
%% BP+SOMF+FK method 
% inds=20;
trace2 = d_bpsomffk(:,inds);
%figure;
% subplot(5,1,3);
plot(t,trace2,'g','linewidth',5); hold on
legend('BP+SOMF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-90 90])
% text(-0.1,125,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',30,'Fontweight','bold');
%% AFK
load d_afk2.mat 

trace3 = d_afk2(:,inds);
% subplot(5,1,4);
plot(t,trace3,'yellow','linewidth',5); hold on
legend('AFK');
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-90 90])
% text(-0.1,125,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',30,'Fontweight','bold');
%% Proposed LO method
% inds=30;
trace4 = d_LO(:,inds);
% subplot(5,1,5)
plot(t,trace4,'r','linewidth',5); hold on
legend('LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-50 100])
% text(-0.1,125,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',30,'Fontweight','bold');
legend('Raw','BP+MF+FK','BP+SOMF+FK','AFK','LO');
annotation(gcf,'line',[0.2609375 0.26015625],...
    [0.178071539657854 0.643856920684292],...
    'Color',[1 0.0745098039215686 0.650980392156863],...
    'LineWidth',8);
annotation(gcf,'textbox',...
    [0.473437500000001 0.577760497667184 0.407812500000001 0.0544323483670296],...
    'String','                              P wave arrival',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'line',[0.640625000000001 0.684375000000002],...
    [0.614307931570762 0.614307931570762],...
    'Color',[1 0.0745098039215686 0.650980392156863],...
    'LineWidth',8);
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
d=20*log(rms(trace4))/log(rms(trace))
end
%% Plot figures

t=[0:n1]*0.0005;
ngap=50;
x=1:n2*3+2*ngap;
% d_1=[din,zeros(n1,ngap0),d1,din-d1,zeros(n1,ngap),d2,din-d2];
% d_2=[d3,din-d3,zeros(n1,ngap),d4,din-d4,zeros(n1,ngap),d5,din-d5];
%% plot the waveforms

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,3,1);LO_imagesc(eq,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

subplot(2,3,2);LO_imagesc(d1_bpmffk,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.518604651162791 0.482945736434109],...
    [0.859808709175739 0.888802488335925],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'textbox',...
    [0.445186046511628 0.587091757387247 0.143961240310077 0.0629860031104205],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

subplot(2,3,3);LO_imagesc(d_bpsomffk,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.794573643410853 0.75891472868217],...
    [0.857475894245723 0.88646967340591],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'textbox',...
    [0.725806201550387 0.589424572317263 0.143961240310078 0.0629860031104205],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

load d_afk2.mat 

subplot(2,3,4);LO_imagesc(d_afk2,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'AFK','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.235658914728682 0.2],...
    [0.383914463452566 0.412908242612753],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15);
annotation(gcf,'textbox',...
    [0.169169350029815 0.294144433892633 0.143961240310077 0.0629860031104204],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

subplot(2,3,5);LO_imagesc(d4,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.522480620155038 0.486821705426356],...
    [0.387024883359253 0.41601866251944],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Improved','signal','energy'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',15);
%% Local similarity
% 
rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(eq-d_dct,d_dct,rect,niter,eps,verb);
[simi2]=LO_localsimi(eq-d_afk2,d_afk2,rect,niter,eps,verb);
[simi3]=LO_localsimi(eq-d1_bpmffk,d1_bpmffk,rect,niter,eps,verb);
[simi4]=LO_localsimi(eq-d_bpsomffk,d_bpsomffk,rect,niter,eps,verb);
[simi5]=LO_localsimi(eq-d4,d4,rect,niter,eps,verb);
% % %
% d_sim1=[simi1,zeros(n1,ngap),simi2,zeros(n1,ngap),simi3];
% d_sim2=[simi4,zeros(n1,ngap),simi5];
%% 
inds1=800:1000;

% dn1=eq-d_dct;
dn2=eq-d_afk2;
dn3=eq-d1_bpmffk;
dn4=eq-d_bpsomffk;
dn5=eq-d4;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
ax(1)= subplot(2,4,1);LO_imagesc(dn3,100,2,x,t);colormap(ax(1),LO_seis);caxis([-25,25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

% annotation(gcf,'textbox',...
ax(2)= subplot(2,4,2);LO_imagesc(dn4,100,2,x,t);colormap(ax(1),LO_seis);caxis([-25,25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax(3)= subplot(2,4,3);LO_imagesc(dn2,100,2,x,t);colormap(ax(1),LO_seis); caxis([-25,25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'AFK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

ax(4)= subplot(2,4,4);LO_imagesc(dn5,100,2,x,t);colormap(ax(1),LO_seis); caxis([-25,25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% ax(5)= subplot(2,4,5);LO_imagesc(simi3,100,1,x,t);colormap(ax(5),jet);caxis([0,0.25]); hold on
ax(1)=subplot(2,4,5); imagesc(x,t,simi3); hold on;
colormap(ax(1),jet); colorbar; caxis([0,0.25]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Channel','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(e)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% colorbar; caxis([0,0.25]);
% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% ax(6)= subplot(2,4,6);LO_imagesc(simi4,100,1,x,t);colormap(ax(6),jet); caxis([0,0.25]); hold on
ax(2)=subplot(2,4,6); imagesc(x,t,simi4); hold on;
colormap(ax(2),jet); colorbar; caxis([0,0.25]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Channel','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(f)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colorbar; caxis([0,0.25]);

% ax(7)= subplot(2,4,7);LO_imagesc(simi2,100,1,x,t);colormap(ax(7),jet);caxis([0,0.25]); hold on
% colormap(ax(7),jet);
ax(3)=subplot(2,4,7); imagesc(x,t,simi2); hold on;
colormap(ax(3),jet); colorbar; caxis([0,0.25]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Channel','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');
text(n2/0.65,-0.04,'AFK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-400,-0.05,'(g)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
colorbar; caxis([0,0.25]);

% ax(8)= subplot(2,4,8);LO_imagesc(simi5,100,1,x,t);colormap(ax(8),jet); caxis([0,0.25]); hold on
ax(4)=subplot(2,4,8); imagesc(x,t,simi5); hold on;
colormap(ax(4),jet); colorbar; caxis([0,0.25]);
ylabel('Time (s)','Fontsize',8,'fontweight','bold');
xlabel('Channel','Fontsize',8,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',8,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(h)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% colorbar; caxis([0, 0.25]);


% % subplot(2,5,5);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));caxis([-25,25]); 
% % ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% % xlabel('Channel','Fontsize',10,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% % text(n2/0.6,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,0.393,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% % subplot(2,6,6);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));caxis([-25,25]); 
% % ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% % xlabel('Channel','Fontsize',10,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% % text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,0.393,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% %%

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,5,6);LO_imagesc(simi1(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,0.397,'DCT','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,4,5);LO_imagesc(simi2(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% subplot(2,4,5);LO_imagesc(simi2,100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,0.397,'SGK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,4,6);LO_imagesc(simi3,100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,0.397,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(2,4,7);LO_imagesc(simi4,100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,0.397,'BP+SOMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,4,8);LO_imagesc(simi5,100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',12,'fontweight','bold');
% xlabel('Channel','Fontsize',12,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
% text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(h)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% 
% % % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% % figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% % subplot(3,2,1);LO_imagesc(simi1(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% % xlabel('Channel','Fontsize',10,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% % text(n2/0.7,-0.04,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% % caxis([0,0.25]);
% % 
% % subplot(3,2,2);LO_imagesc(simi2(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% % xlabel('Channel','Fontsize',14,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% % text(n2/0.7,-0.04,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% % caxis([0,0.25]);
% % 
% % subplot(3,2,3);LO_imagesc(simi3(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% % xlabel('Channel','Fontsize',14,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% % text(n2/0.7,-0.04,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% % caxis([0,0.25]);
% % 
% % % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% % subplot(3,2,4);LO_imagesc(simi4(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% % xlabel('Channel','Fontsize',14,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% % text(n2/0.7,-0.04,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% % caxis([0,0.25]);
% % 
% % subplot(3,2,5);LO_imagesc(simi5(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% % ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% % xlabel('Channel','Fontsize',14,'fontweight','bold');
% % set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% % text(n2/0.65,-0.04,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% % text(-300,-0.1,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% % c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% % caxis([0,0.25]);
% % 
% % print(gcf,'-depsc','-r200','fig10.eps');
% % 
% % 
% % 
% % 
% a =
% 
%    13.3315
% 
% 
% b =
% 
%    14.4772
% 
% 
% c =
% 
%    15.1647
%% 
% n1 = 2000; n2 = 960
%4.353120 seconds.
%489.525771 seconds.
% 1113.018505 seconds.

% 1000 960 
% 3.89
% 497.76
% 863.63
% 2000 and 1500
