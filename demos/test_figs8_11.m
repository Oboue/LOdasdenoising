% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example
%% first real example 
% Script to plot Figures 1 and 2

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
% addpath(genpath('./MATseisdl'));
% addpath(genpath('./subroutines'));
%% addpath
addpath('../MATseisdl/');
addpath('../subroutines/');
addpath('../LO_data/');
addpath('../LOsrcs/');
addpath('../seistr/');
%% Load the raw DAS seismic data corrupted by strong hogh-frequency noise, high-amplitude 
%% Load the DAS data

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
ngap=50;
x=1:n2*3+2*ngap;

for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('LO_data/eq-',num2str(ii),'.mat')
        load(strcat('LO_data/eq-',num2str(ii),'.mat'));
    end
%   d1=d1;
    eq=d1;
    din=d1;

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
d_bpsomffk=d1-LO_fk_dip(d1,0.02);%
toc
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
par.ns=8;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % size of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.08;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=.15;              % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 5;
par.rec(2) = 5;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;                 % verbosity flag (default: 0) 

%
tic
d_LO=LO(eq,par);
d4=d_LO;
toc
%
end
%% Plot figures
% eq=zeros(2000,960);
% [n1,n2]=size(eq);
%% Raw data
t=[0:n1-1]*0.0005;
inds=30;
trace = eq(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); hold on
subplot(4,1,1);
plot(t,trace,'black','linewidth',6); hold on
legend('Raw');
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-300 500])
text(-0.08,600,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% BP+SOMF+FK method
inds=30;
trace1 = d1_bpmffk(:,inds);
%figure;
subplot(4,1,2);
plot(t,trace1,'b','linewidth',6); hold on
legend('BP+MF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-300 500])
 text(-0.08,600,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% BP+SOMF+FK method 
inds=30;
trace2 = d_bpsomffk(:,inds);
%figure;
subplot(4,1,3);
plot(t,trace2,'g','linewidth',6); hold on
legend('BP+SOMF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-300 500])
text(-0.08,600,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% Proposed LO method
inds=30;
trace3 = d_LO(:,inds);
subplot(4,1,4)
plot(t,trace3,'r','linewidth',6); hold on
legend('LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-300 500])
text(-0.08,600,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
%% plot the waveforms
t=[0:n1]*0.0005;
% length(t)
ngap=50;
x=1:n2*3+2*ngap;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,2,1);LO_imagesc(eq,100,2,x,t); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.132222222222222 0.754679802955665 0.332222222222222 0.0328339412257772],...
    'Color',[1 0 0],...
    'LineWidth',5,...
    'FontWeight','bold',...
    'FontSize',3,...
    'FitBoxToText','off',...
    'EdgeColor',[1 0 0]);

subplot(2,2,2);LO_imagesc(d1_bpmffk,100,2,x,t); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.573333333333333 0.756650246305416 0.332222222222222 0.0328339412257772],...
    'Color',[1 0 0],...
    'LineWidth',5,...
    'FontWeight','bold',...
    'FontSize',3,...
    'FitBoxToText','off',...
    'EdgeColor',[1 0 0]);

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,2,3);LO_imagesc(d_bpsomffk,100,2,x,t); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.131111111111111 0.279802955665022 0.332222222222222 0.0328339412257773],...
    'Color',[1 0 0],...
    'LineWidth',5,...
    'FontWeight','bold',...
    'FontSize',3,...
    'FitBoxToText','off',...
    'EdgeColor',[1 0 0]);

subplot(2,2,4);LO_imagesc(d4,100,2,x,t); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.571111111111111 0.279802955665023 0.332222222222222 0.0328339412257772],...
    'Color',[1 0 0],...
    'LineWidth',5,...
    'FontWeight','bold',...
    'FontSize',3,...
    'FitBoxToText','off',...
    'EdgeColor',[1 0 0]);
%% Zoomed sections 

inds1=800:1000;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,2,1);LO_imagesc(eq(inds1,:),100,2,x,t(inds1)); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.380555555555556 0.399444444444444],...
    [0.781590789097661 0.742259422332667],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',20);

subplot(2,2,2);LO_imagesc(d1_bpmffk(inds1,:),100,2,x,t(inds1)); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.821666666666667 0.840555555555555],...
    [0.777245932448864 0.737914565683869],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',20);

subplot(2,2,3);LO_imagesc(d_bpsomffk(inds1,:),100,2,x,t(inds1)); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.382777777777778 0.401666666666666],...
    [0.30885972942733 0.269528362662336],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',20);

subplot(2,2,4);LO_imagesc(d4(inds1,:),100,2,x,t(inds1)); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.820555555555556 0.839444444444444],...
    [0.30729131852111 0.267959951756116],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',20,...
    'HeadLength',20);
%% local similarity metrics
% 
rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(eq-d_dct,d_dct,rect,niter,eps,verb);
% [simi2]=LO_localsimi(eq-d11,d11,rect,niter,eps,verb);
[simi3]=LO_localsimi(eq-d1_bpmffk,d1_bpmffk,rect,niter,eps,verb);
[simi4]=LO_localsimi(eq-d_bpsomffk,d_bpsomffk,rect,niter,eps,verb);
[simi5]=LO_localsimi(eq-d4,d4,rect,niter,eps,verb);
%
% d_sim1=[simi1,zeros(n1,ngap),simi2,zeros(n1,ngap),simi3];
% d_sim2=[simi4,zeros(n1,ngap),simi5];
% % 
inds1=800:1000;

% dn1=eq-d_dct;
% dn2=eq-d11;
dn3=eq-d1_bpmffk;
dn4=eq-d_bpsomffk;
dn5=eq-d4;

figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');

subplot(2,3,1);LO_imagesc(dn3(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

subplot(2,3,2);LO_imagesc(dn4(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+SOSVMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

subplot(2,3,3);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
 
subplot(2,3,4);LO_imagesc(simi3(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,3,5);LO_imagesc(simi4(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');

subplot(2,3,6);LO_imagesc(simi5(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-300,0.393,'(f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');


% % % 
% % % 
% % % 
% % % SNR =
% a =
% 
%    18.3820
% 
% 
% b =
% 
%    18.9864
% 
% 
% c =
% 
%    19.5895

% size n1 = 2000; n2 = 960 
% BP+MF+FK= 4.089881 seconds.
% BP+SOMF+FK = 479.896468 seconds
% proposed 877.799570 seconds.