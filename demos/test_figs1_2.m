% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

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
addpath('../LO_data/');
addpath('../LOsrcs/');
addpath('../seistr/');
%% Load the raw DAS seismic data corrupted by strong hogh-frequency noise, high-amplitude 

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

for ii=4
    
if ~ismember(NOs(ii),[14,16,17,27,47,52])
    load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
end
eq=d1;
din=d1;
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
tic
d1=LO_bandpass(din,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
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
ns=15;                         % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
tic
d2=LO_bandpasssosvmf(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
%% Denoising using the BP+SOSVMF+FK method 
%  Parameter tuning：add the key parameters of the dip filter in FK domain method
%
w=0.08;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.5;             % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
% 
tic
d4=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
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
rect(2)=rect(2);            % "      "        "
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
par.niter1=niter1;           % Number of iteration
%
rec = zeros(3, 1);          % 3-D vector denoting smooth radius 
par.rec(1) = 150;
par.rec(2) = 150;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb=1;                 % verbosity flag (default: 0) 

%
tic
d_LO=LO(din,par);
d5=d_LO;
toc
% 
end
%% Plot figures

t=[0:n1]*0.0005;
ngap0=1000;

ngap=50;
x=1:n2*3+2*ngap;
d_1=[din,zeros(n1,ngap0),d1,din-d1,zeros(n1,ngap),d2,din-d2];
d_2=[d3,din-d3,zeros(n1,ngap),d4,din-d4,zeros(n1,ngap),d5,din-d5];

% combined figure
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);LO_imagesc(d_1,100,2,x,t);caxis([-25,25]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/3.5,-0.03,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.79,-0.03,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,-0.03,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.12995951417004 0.754408060453399 0.7758987854251 0.0340050377833753],...
    'Color',[1 0 0],...
    'LineWidth',3);
annotation(gcf,'textarrow',[0.468608293655442 0.497354497354497],...
    [0.679535637681291 0.635036496350365],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.708026283073435 0.728835978835977],...
    [0.672265823903992 0.614384687439247],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Honrizontal','noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

subplot(2,1,2);LO_imagesc(d_2,100,2,x,t);caxis([-20,20]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/2.0,-0.03,'BP+SOSVMF+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.64,-0.03,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/0.65,-0.03,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.135870758987019 0.14374695863747 0.116774743658483 0.0722270315960674],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.459426347338605 0.427248677248677],...
    [0.210804497048998 0.237575757575758],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Weaker','signal energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.713262326174592 0.69179894179894],...
    [0.206851587471678 0.243636363636364],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'rectangle',...
    [0.129959514170039 0.280468666514005 0.775898785425098 0.0340050377833752],...
    'Color',[1 0 0],...
    'LineWidth',3);
print(gcf,'-depsc','-r300','fig1.eps');
%% Zoomed sections 

inds1=800:1000;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);LO_imagesc(d_1(inds1,:),100,2,x,t(inds1));caxis([-20,20]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
annotation(gcf,'textarrow',[0.468608293655442 0.497354497354497],...
    [0.679535637681291 0.635036496350365],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.708026283073435 0.728835978835977],...
    [0.678319092669125 0.62043795620438],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Honrizontal','noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
subplot(2,1,2);LO_imagesc(d_2(inds1,:),100,2,x,t(inds1));caxis([-20,20]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
annotation(gcf,'textbox',...
    [0.13719351030977 0.197080291970803 0.116774743658483 0.0722270315960675],...
    'Color',[1 0 0],...
    'String',{'Random','background','noise'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textarrow',[0.455458093370352 0.485449735449735],...
    [0.208380254624756 0.173965936739659],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Weaker','signal','energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.715907828820097 0.746031746031746],...
    [0.214124314744405 0.171532846715329],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal','energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
text(-200,0.25,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3.5,0.255,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.65,0.255,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,0.255,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-200,0.395,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.0,0.395,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.65,0.395,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.38,0.395,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');

print(gcf,'-depsc','-r300','fig2.eps');

%%