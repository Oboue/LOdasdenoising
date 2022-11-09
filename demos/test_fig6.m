% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figure 6

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
%% Ccomparison between SOMF, SOSVMF methods
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
[n1,n2]=size(eq);
t=[0:n1]*0.0005;
ngap=50;
x=1:n2*3+2*ngap;

for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('LO_data/eq-',num2str(ii),'.mat')
        load(strcat('LO_data/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
    din=d1;
%% Denoising using the BP+MF+FK+Curvelet method

dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
ns=8;      % spray radius
nfw=ns*2+1;
ifb=1;
axis=2;
%
w=0.05;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;                % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
%
tic
d1=LO_bpmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,nfw,ifb,axis,w,n1,n2,c1,c2,c3,niter1);
toc

%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 30;
rec(2) = 30;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb2=1;              % verbosity flag (default: 0) 

%
tic
d2=LO_bpmffkcurveletlow(din,dt,flo,fhi,nplo,nphi,phase,verb0,nfw,ifb,axis,w,n1,n2,c1,c2,c3,niter1,rec,eps1,niter2,verb2);
toc
%

comp1=[eq,zeros(n1,ngap),d1,eq-d1,zeros(n1,ngap),d2,eq-d2]; 
%% LO method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag
%
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
w=0.05;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;                % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

d3=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
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
par.rec(1) = 30;
par.rec(2) = 30;
par.rec(3) = 1;
par.eps1=0;                 % regularization parameter, default 0.0
par.niter2=20;              % number of CG iterations
par.verb2=1;                % verbosity flag (default: 0) 
%
tic
d4=LO(din,par);
toc
%
comp2=[eq,zeros(n1,ngap),d3,eq-d3,zeros(n1,ngap),d4,eq-d4]; 
end
%%
%% Plot figures 

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);LO_imagesc(comp1,95,1,x,t);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
text(n2/3,-0.04,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/4+ngap+n2,-0.04,'BP+MF+FK+Curvelet','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5+ngap*2+n2*2,-0.04,'BP+MF+FK+Curvelet+LOW (LO)','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.791005291005289 0.775132275132273],...
    [0.789537712895378 0.824817518248176],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.776455026455025 0.772486772486769],...
    [0.892944038929441 0.862530413625305],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.488095238095235 0.465608465608462],...
    [0.771289537712897 0.811435523114356],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.460317460317458 0.456349206349202],...
    [0.884428223844283 0.854014598540147],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

subplot(2,1,2);LO_imagesc(comp2,95,1,x,t);
ylabel('Time (s)','Fontsize',13,'fontweight','bold');
xlabel('Channel','Fontsize',13,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',13,'Fontweight','bold');
text(n2/3,-0.04,'Raw data','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/4+ngap+n2,-0.04,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5+ngap*2+n2*2,-0.04,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.48148148148148 0.465608465608464],...
    [0.327250608272506 0.362530413625304],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Weaker','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.791005291005289 0.775132275132273],...
    [0.324817518248175 0.360097323600973],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

print(gcf,'-depsc','-r300','fig6.eps');


