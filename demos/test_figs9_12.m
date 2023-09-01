% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figures 9, 10, 11 and 12

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
% subplot(5,1,1);
plot(t,trace,'black','linewidth',6); hold on
legend('Raw');
xlabel('Time (s)')
ylabel('Amplitude')
ylim([-300 500])
text(-0.08,600,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',40,'Fontweight','bold');
%% BP+SOMF+FK method
inds=30;
trace1 = d1_bpmffk(:,inds);
%figure;
% subplot(5,1,2);
plot(t,trace1,'b','linewidth',6); hold on
legend('BP+MF+FK');
xlabel('Time (s)');
ylabel('Amplitude');
 ylim([-300 500])
 text(-0.08,600,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',40,'Fontweight','bold');
%% BP+SOMF+FK method 
inds=30;
trace2 = d_bpsomffk(:,inds);
%figure;
% subplot(5,1,3);
plot(t,trace2,'g','linewidth',6); hold on
legend('BP+SOMF+FK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-200 500])
text(-0.08,600,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',40,'Fontweight','bold');
%% AFK
load d_afk1 % load afk data
inds=30;
trace3 = d_afk1(:,inds);
%figure;
% subplot(5,1,4);
plot(t,trace3,'yellow','linewidth',6); hold on
legend('AFK');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-300 500])
text(-0.08,600,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',40,'Fontweight','bold');
%% Proposed LO method
inds=30;
trace4 = d_LO(:,inds);
% subplot(5,1,5)
plot(t,trace4,'r','linewidth',6); hold on
legend('LO');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-300 500])
legend('Raw','BP+MF+FK','BP+SOMF+FK','AFK','LO');
text(-0.08,600,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',4,'Fontsize',30,'Fontweight','bold');
annotation(gcf,'textbox',...
    [0.473437500000001 0.577760497667184 0.407812500000001 0.0544323483670296],...
    'String','                              P wave arrival',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'line',[0.640625000000001 0.684375000000002],...
    [0.614307931570762 0.614307931570762],'Color',[0 1 1],'LineWidth',8);
annotation(gcf,'line',[0.323437500000001 0.322656250000001],...
    [0.209175738724728 0.674961119751166],'Color',[0 1 1],'LineWidth',8);
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
d=20*log(rms(trace4))/log(rms(trace))
%% plot the waveforms
t=[0:n1]*0.0005;
% length(t)
ngap=50;
x=1:n2*3+2*ngap;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,3,1);LO_imagesc(eq,100,2,x,t);caxis([-60,60]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'Raw','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.12990625 0.755832037325039 0.21228125 0.0318818040435459],...
    'Color',[1 0 0],...
    'LineWidth',5);

subplot(2,3,2);LO_imagesc(d1_bpmffk,100,2,x,t);caxis([-60,60]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+MF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.4119375 0.757387247278383 0.21071875 0.031881804043546],'Color',[1 0 0],...
    'LineWidth',5);

% 
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,3,3);LO_imagesc(d_bpsomffk,100,2,x,t);caxis([-60,60]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'BP+SOMF+FK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.692406250000002 0.755832037325039 0.21071875 0.031881804043546],...
    'Color',[1 0 0],...
    'LineWidth',5);

subplot(2,3,4);LO_imagesc(d_afk1,100,2,x,t);caxis([-60,60]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,caxis([-40,40]); hold ont(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'AFK','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.130687500000001 0.28149300155521 0.21071875 0.0318818040435461],...
    'Color',[1 0 0],...
    'LineWidth',5);

subplot(2,3,5);LO_imagesc(d4,100,2,x,t);caxis([-60,60]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,-0.04,'LO','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
text(-250,-0.05,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.4119375 0.280715396578538 0.21071875 0.031881804043546],'Color',[1 0 0],...
    'LineWidth',5);
%% Zoomed sections 

inds1=800:1000;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,3,1);LO_imagesc(eq(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'Raw','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.25859375 0.3015625],...
    [0.774272161741835 0.746500777604977],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',25);

subplot(2,3,2);LO_imagesc(d1_bpmffk(inds1,:),100,2,x,t(inds1)); caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'BP+MF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.541406250000001 0.584375000000001],...
    [0.775049766718507 0.747278382581649],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',25);

subplot(2,3,3);LO_imagesc(d_bpsomffk(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'BP+SOMF+FK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.821093750000001 0.864062500000001],...
    [0.771161741835148 0.74339035769829],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',25);

subplot(2,3,4);LO_imagesc(d_afk1(inds1,:),100,2,x,t(inds1)); caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'AFK','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.260156250000001 0.303125000000001],...
    [0.303821150855365 0.276049766718507],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',25);

subplot(2,3,5);LO_imagesc(d4(inds1,:),100,2,x,t(inds1)); caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/0.65,0.396,'LO','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.540625000000001 0.583593750000001],...
    [0.292157076205288 0.264385692068429],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',25);
%% local similarity metrics
% 
rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(eq-d_dct,d_dct,rect,niter,eps,verb);
% [simi2]=LO_localsimi(eq-d11,d11,rect,niter,eps,verb);
[simi3]=LO_localsimi(eq-d1_bpmffk,d1_bpmffk,rect,niter,eps,verb);
[simi4]=LO_localsimi(eq-d_bpsomffk,d_bpsomffk,rect,niter,eps,verb);
[simi5]=LO_localsimi(eq-d_afk1,d_afk1,rect,niter,eps,verb);
[simi6]=LO_localsimi(eq-d4,d4,rect,niter,eps,verb);
%
% d_sim1=[simi1,zeros(n1,ngap),simi2,zeros(n1,ngap),simi3];
% d_sim2=[simi4,zeros(n1,ngap),simi5];
% % 
inds1=800:1000;

% dn1=eq-d_dct;
% dn2=eq-d11;
dn3=eq-d1_bpmffk;
dn4=eq-d_bpsomffk;
dn5=eq-d_afk1;
dn6=eq-d4;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,4,1);LO_imagesc(dn3(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,0.397,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(a)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

subplot(2,4,2);LO_imagesc(dn4(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,0.397,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(b)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

subplot(2,4,3);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,0.397,'AFK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(c)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
 
subplot(2,4,4);LO_imagesc(dn6(inds1,:),100,2,x,t(inds1));caxis([-60,60]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,0.397,'LO','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-500,0.393,'(d)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
 
% subplot(2,4,5);LO_imagesc(simi3(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ax(1)=subplot(2,4,5); imagesc(x,t,simi3(inds1,:)); hold on;
colormap(ax(1),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.03,'BP+MF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-1000,-0.04,'(e)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(2,4,6);LO_imagesc(simi4(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ax(2)=subplot(2,4,6); imagesc(x,t,simi4(inds1,:)); hold on;
colormap(ax(2),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.03,'BP+SOMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-1000,-0.04,'(f)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(2,4,7);LO_imagesc(simi5(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ax(3)=subplot(2,4,7); imagesc(x,t,simi5(inds1,:)); hold on;
colormap(ax(3),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.03,'AFK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-1000,-0.04,'(g)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

% subplot(2,4,8);LO_imagesc(simi6(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.5]); hold on
ax(4)=subplot(2,4,8); imagesc(x,t,simi6(inds1,:)); hold on;
colormap(ax(4),jet); colorbar; caxis([0,0.5]);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.03,'LO','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-1000,-0.04,'(h)','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');

