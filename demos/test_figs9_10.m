% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figures 9 and 10

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
%   d1=d1;
    eq=d1;
    din=d1;

%% Denosing using the SOSVMF+FK+Curvelet method
%  Parameter tuning for the BP method

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

w=0.05;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;                % Thresholding parameter (alpha)
niter1=10;           % Number of iteration
%
tic
d_sosvmffkcurvelet=LO_sosvmffkcurvelet(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%

% Denosing using the SOSVMF+FK+Curvelet+LOW (LO) method

rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 30;
rec(2) = 30;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb2=1;              % verbosity flag (default: 0)
%
tic
d_sosvmffkcurveletlow=LO_sosvmffkcurveletlow(din,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1,rec,eps1,niter2,verb2);
toc
%
comp1=[din,zeros(n1,ngap),d_sosvmffkcurvelet,din-d_sosvmffkcurvelet,zeros(n1,ngap),d_sosvmffkcurveletlow,din-d_sosvmffkcurveletlow]; 
%% BP+FK+Curvelet method
%
dt=0.0005; % sampling
flo=0;     % Low frequency in band, default is 0
fhi=200;   % High frequency in band, default is Nyquist
nplo=6;    % number of poles for low cutoff
nphi=6;    % number of poles for high cutoff
phase=0;   % y: minimum phase, n: zero phase
verb0=0;   % verbosity flag

w=0.02;              % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=2.5;                % Thresholding parameter (alpha)
niter1=10;           % Number of iteration

tic
d_bpfkcurvelet=LO_bpfkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,w,n1,n2,c1,c2,c3,niter1);
toc
%% BP+FK+Curvelet+LOW (LO) method

rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
rec(1) = 30;
rec(2) = 30;
rec(3) = 1;
eps1=0;               % regularization parameter, default 0.0
niter2=20;            % number of CG iterations
verb2=1;              % verbosity flag (default: 0)

d_bpfkcurveletlow=LO_bpfkcurveletlow(din,dt,flo,fhi,nplo,nphi,phase,verb0,w,n1,n2,c1,c2,c3,niter1,rec,eps1,niter2,verb2);
%%
comp2=[din,zeros(n1,ngap),d_bpfkcurvelet,din-d_bpfkcurvelet,zeros(n1,ngap),d_bpfkcurveletlow,din-d_bpfkcurveletlow]; 
%% BP+SOSVMF+FK method
tic
d_bpsosvmffk=LO_bpsosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% BP+SOSVMF+FK+LOW (LO) method
tic
d_bpsosvmffklow=LO_bpsosvmffklow(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,rec,niter2,eps1,verb2);
toc
%%
comp3=[din,zeros(n1,ngap),d_bpsosvmffk,din-d_bpsosvmffk,zeros(n1,ngap),d_bpsosvmffklow,din-d_bpsosvmffklow]; 

%% BP+SOSVMF+FK+Curvelet method
tic
d_bandpasssosvmffkcurvelet=LO_bandpasssosvmffkcurvelet(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w,c1,c2,c3,niter1);
toc
%
%% BP+SOSVMF+FK+Curvelet+LOW (LO) method
%% LO method
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
par.w=0.02;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=2.5;               % Thresholding parameter (alpha)
par.niter1=10;           % Number of iteration
%
rec = zeros(3, 1);    % 3-D vector denoting smooth radius 
par.rec(1) = 30;
par.rec(2) = 30;
par.rec(3) = 1;
par.eps1=0;               % regularization parameter, default 0.0
par.niter2=20;            % number of CG iterations
par.verb=1;               % verbosity flag (default: 0) 
%
tic
d_bandpasssosvmffkcurveletlow=LO(eq,par);
toc
%
%%
comp4=[din,zeros(n1,ngap),d_bandpasssosvmffkcurvelet,din-d_bandpasssosvmffkcurvelet,zeros(n1,ngap),d_bandpasssosvmffkcurveletlow,din-d_bandpasssosvmffkcurveletlow]; 
end
%%
t=[0:n1]*0.0005;
% ngap0=1000;

ngap=50;
x=1:n2*3+2*ngap;

% combined figure
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(4,1,1);LO_imagesc(comp1,100,2,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/3.5,-0.09,'Raw data','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.79,-0.09,'SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,-0.09,'SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.361111111111111 0.407407407407407],...
    [0.861364087765208 0.836538461538462],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.518518518518516 0.558201058201058],...
    [0.859236260398166 0.836538461538462],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal','leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.685185185185184 0.723544973544973],...
    [0.86022788344705 0.836538461538463],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.830687830687831 0.87037037037037],...
    [0.852163461538462 0.835336538461539],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
subplot(4,1,2);LO_imagesc(comp2,100,2,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/3.5,-0.09,'Raw data','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.79,-0.09,'BP+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,-0.09,'BP+FK+Curvelet+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(b)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.363756613756613 0.411375661375661],...
    [0.641412164688286 0.671875],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.535714285714285 0.563492063492063],...
    [0.667067307692308 0.693509615384617],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.69973544973545 0.722222222222222],...
    [0.634615384615385 0.661057692307693],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.846560846560847 0.873015873015871],...
    [0.661057692307692 0.688701923076924],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

subplot(4,1,3);LO_imagesc(comp3,100,2,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/3.5,-0.09,'Raw data','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.79,-0.09,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,-0.09,'BP+SOSVMF+FK+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(c)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.477513227513225 0.461640211640211],...
    [0.420673076923078 0.451923076923078],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Visible','signal','leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.791005291005288 0.775132275132274],...
    [0.426682692307693 0.457932692307693],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

subplot(4,1,4);LO_imagesc(comp4,100,2,x,t);
ylabel('Time (s)','Fontsize',12,'fontweight','bold');
xlabel('Channel','Fontsize',12,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',12,'Fontweight','bold');
text(n2/3.5,-0.09,'Raw data','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.79,-0.09,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,-0.09,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',...
    [0.12995951417004 0.20818957305844 0.7758987854251 0.0340050377833751],...
    'Color',[1 0 0],...
    'LineWidth',3);

print(gcf,'-depsc','-r200','fig9.eps');

%%  Zoom-in section
inds1=400:800;

%local similarity maps

rect=[30,30,1];niter=20;eps=0;verb=0;
[simi1]=LO_localsimi(din-d_bandpasssosvmffkcurvelet,d_bandpasssosvmffkcurvelet,rect,niter,eps,verb);
[simi2]=LO_localsimi(din-d_bandpasssosvmffkcurveletlow,d_bandpasssosvmffkcurveletlow,rect,niter,eps,verb);
%%
d_sim=[simi1,zeros(n1,ngap),simi2];

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);LO_imagesc(comp4(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
annotation(gcf,'textarrow',[0.390713964398158 0.412698412698413],...
    [0.683603092317042 0.646844660194174],'Color',[1 0 0],...
    'String',{'Weaker','   signal energy'},...
    'LineWidth',2,...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.482160517144917 0.465555555555556],...
    [0.8123458831972 0.837162837162837],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.79201058201058 0.777460317460317],...
    [0.815881933600379 0.846221739425622],'Color',[1 0 0],...
    'String',{'Weaker','signal leakage'},...
    'LineWidth',2,...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.688190335280036 0.723544973544973],...
    [0.675480725700934 0.643203883495145],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal energy'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
text(-200,0.19,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/3,0.19,'Raw data','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.8,0.19,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.4,0.19,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');

print(gcf,'-depsc','-r300','fig10.eps');
%%
%figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,2);LO_imagesc(d_sim(inds1,:),100,2,x,t(inds1));colormap(jet);hold on
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(-200,0.19,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.5,0.19,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(n2/0.43,0.19,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',13,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.246666666666667 0.196666666666667],...
    [0.292878635907723 0.373119358074223],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FontName','Helvetica Neue');
annotation(gcf,'arrow',[0.23 0.15],...
    [0.223671013039117 0.164493480441324],'Color',[1 0 0],'LineWidth',2);

c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 13;%c.Label.FontWeight = bold;
caxis([0,1]);

print(gcf,'-depsc','-r200','fig10.eps');