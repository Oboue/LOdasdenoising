% Protecting the weak signals in distributed acoustic sensing data processing using local orthogonalization: the FORGE data example

% Script to plot Figure 4

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
%----------------------------------------------------------------------------------------------------------
clc;clear;close all;
%% addpath
addpath('../LO_data/');
addpath('../LOsrcs/');
addpath('../seistr/');
%% Ccomparison between MF, SOMF, and SOSVMF methods
%% Load the DAS data

eq=zeros(2000,960);
[n1,n2]=size(eq);
for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('LO_data/eq-',num2str(ii),'.mat')
        load(strcat('LO_data/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
%% Denosing using the MF and MF+LOW methods 
ns=8;                         % spray radius
%
tic
d_mf=LO_mf(eq,ns*2+1,1,2);
toc
%
noi1=eq-d_mf;
% prepare paramters for ortho
rect = zeros(3,1);
rect(1) = 30;
rect(2) = 30;
rect(3) = 1;
eps=0;
niter=20;
verb=1; 

[d_mflow,noi2,low]=LO_localortho(d_mf,noi1,rect,niter,eps,verb);  
%% Denoising using the SOMF and SOMF+LOW methods

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=100;                   % rect:  smoothing radius (ndim*1)
rect(2)=100;                   % "      "        "
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
d_somf=LO_somf(eq,niter,liter,order1,eps_dv, eps_cg, tol_cg,rect,verb,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc

noi1=eq-d_somf;
% prepare paramters for ortho
rect = zeros(3,1);
rect(1) = 30;
rect(2) = 30;
rect(3) = 1;
eps=0;
niter=20;
verb=1;

[d_somflow,noi2,low]=LO_localortho(d_somf,noi1,rect,niter,eps,verb);  
%
%% Denosing using the SOSVMF and SOSVMF+LOW method 
% Parameter tuning：add the key parameters of the SOSVMF method

niter=2;                      % number of nonlinear iterations
liter=10;                     % liter: number of linear iterations (in divn)
order1=3;                     % order: accuracy order
eps_dv=0.01;                  % eps_dv: eps for divn  (default: 0.01)
eps_cg=1;                     % eps_cg: eps for CG    (default: 1)
tol_cg=0.000001;              % tol_cg: tolerence for CG (default: 0.000001)
rect(1)=100;                   % rect:  smoothing radius (ndim*1)
rect(2)=100;                   % "      "        "
rect(3)=1;                    % "      "        "
verb=1;                       % verbosity flag

adj=0;                        % adjoint flag
add=0;                        % adding flag
ns=7;                        % spray radius
order2=2;                     % PWD order
eps=0.01;                     % regularization (default:0.01);
ndn=n1*n2;                    % size of dn (n1*n2)
nds=n1*n2;                    % size of ds (n1*n2)
type_mf=1;                    % 0 (MF) or 1 (SVMF)
ifsmooth=0;                   % 1 (if smooth) or 0 (only MF)
%
tic
d_sosvmf=LO_sosvmf(eq,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth);
toc
%
noi1=eq-d_sosvmf;
% prepare paramters for ortho
rect = zeros(3,1);
rect(1) = 30;
rect(2) = 30;
rect(3) = 1;
eps=0;
niter=20;
verb=1;

[d_sosvmflow,noi2,low]=LO_localortho(d_sosvmf,noi1,rect,niter,eps,verb);  
end
%%

t=[0:n1]*0.0005;
x=1:n2;
ngap=50;

d_1=[eq,zeros(n1,ngap),d_mf,eq-d_mf,zeros(n1,ngap),d_mflow,eq-d_mflow];
d_2=[eq,zeros(n1,ngap),d_somf,eq-d_somf,zeros(n1,ngap),d_somflow,eq-d_somflow];
d_3=[eq,zeros(n1,ngap),d_sosvmf,eq-d_sosvmf,zeros(n1,ngap),d_sosvmflow,eq-d_sosvmflow];
%% figure;
inds1=1:500;

figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(3,1,1);LO_imagesc(d_1(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/9,-0.015,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/200,-0.05,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.489417989417989 0.466931216931217],...
    [0.797330097087378 0.764563106796116],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Stronger','signal','leakage'},...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.808201058201056 0.785714285714284],...
    [0.814320388349514 0.781553398058252],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.808201058201056 0.785714285714284],...
    [0.814320388349514 0.781553398058252],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
% text(0.1,0.95,labels{2},'color','b','Fontsize',13,'fontweight','bold','HorizontalAlignment','left');
subplot(3,1,2);LO_imagesc(d_2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/9,-0.015,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5,-0.015,'SOMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.015,'SOMF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/200,-0.05,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.510582010582008 0.488095238095236],...
    [0.52791262135922 0.495145631067958],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.808201058201058 0.793650793650791],...
    [0.520631067961165 0.490291262135921],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
subplot(3,1,3);LO_imagesc(d_3(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',14,'fontweight','bold');
xlabel('Channel','Fontsize',14,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
text(n2/9,-0.015,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/2.5,-0.015,'SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/1.25,-0.015,'SOSVMF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/200,-0.05,'(c)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.513227513227512 0.490740740740739],...
    [0.231796116504853 0.199029126213591],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');
annotation(gcf,'textarrow',[0.82275132275132 0.800264550264548],...
    [0.239077669902912 0.206310679611649],'Color',[1 0 0],'TextColor',[1 0 0],...
    'LineWidth',2,...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Helvetica Neue');

print(gcf,'-depsc','-r300','fig4.eps');