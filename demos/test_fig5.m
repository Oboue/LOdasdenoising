% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

%  Script to plot Figure 5

%  Copyright (C) Oboue et al., 2022

%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version. b

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
for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('LO_data/eq-',num2str(ii),'.mat')
        load(strcat('LO_data/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
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
toc
%

ns=8;  % spray radius
%
tic
d_mf=LO_mf(d1,ns*2+1,1,2);
toc

%FK
d1_bpmffk=d_mf-LO_fk_dip(d_mf,0.02);%
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
d1_LO=LO(eq,par);
toc
%
end
%%

t=[0:n1]*0.0005;
x=1:n2;
ngap=50;

d_1=[eq,zeros(n1,ngap),d1_bpmffk,eq-d1_bpmffk,zeros(n1,ngap),d1_LO,eq-d1_LO];

for ii=2
    if ~ismember(NOs(ii),[14,16,17,27,47,52])
%         strcat('amf_data/eq-',num2str(ii),'.mat')
        load(strcat('LO_data/eq-',num2str(NOs(ii)),'.mat'));
    end
    d1=d1;
    eq=d1;
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
d2=LO_bandpass(eq,dt,flo,fhi,nplo,nphi,phase,verb0);
toc
%

ns=8;  % spray radius
%
tic
d2_mf=LO_mf(d2,ns*2+1,1,2);
toc

%FK
d2_bpmffk=d2_mf-LO_fk_dip(d2_mf,0.02);%
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
par.ns=10;                         % spray radius
par.order2=2;                     % PWD order
par.eps=0.01;                     % regularization (default:0.01);
par.ndn=n1*n2;                    % size of dn (n1*n2)
par.nds=n1*n2;                    % siz e of ds (n1*n2)
par.type_mf=1;                    % 0 (MF) or 1 (SVMF)
par.ifsmooth=0;                   % 1 (if smooth) or 0 (only MF);
% 
par.w=0.02;   % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
par.c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
par.c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
par.c3=0.25;               % Thresholding parameter (alpha)
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
d2_LO=LO(eq,par);
toc
%
end

d_2=[eq,zeros(n1,ngap),d2_bpmffk,eq-d2_bpmffk,zeros(n1,ngap),d2_LO,eq-d2_LO];
% d_3=[eq,zeros(n1,ngap),d_sosvmf,eq-d_sosvmf,zeros(n1,ngap),d_sosvmflow,eq-d_sosvmflow];
%% figure;
inds1=1:800;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(2,1,1);LO_imagesc(d_1(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/9,-0.015,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5,-0.015,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.015,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/200,-0.05,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.38,labels{6},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'textarrow',[0.50132275132275 0.466931216931216],...
    [0.663824604141291 0.697929354445798],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.810846560846557 0.776455026455023],...
    [0.672350791717417 0.706455542021923],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Weaker','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

inds2=1000:1400;

subplot(2,1,2);LO_imagesc(d_2(inds2,:),100,2,x,t(inds2));caxis([-50,50]);
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/9,0.49,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5,0.49,'BP+MF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,0.49,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/200,0.47,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
text(0.1,0.69,labels{2},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
annotation(gcf,'textarrow',[0.5026455026455 0.568783068783063],...
    [0.352009744214372 0.393422655298416],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.51190476190476 0.563492063492061],...
    [0.236297198538367 0.276492082825821],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.828042328042326 0.878306878306875],...
    [0.361753958587088 0.403166869671132],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.818783068783064 0.870370370370365],...
    [0.248477466504261 0.288672350791715],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Weaker','signal leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

print(gcf,'-depsc','-r300','fig5.eps');


