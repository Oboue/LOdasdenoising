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

% d_1=[eq,zeros(n1,ngap),d_mf,eq-d_mf,zeros(n1,ngap),d_mflow,eq-d_mflow];
% d_2=[eq,zeros(n1,ngap),d_somf,eq-d_somf,zeros(n1,ngap),d_somflow,eq-d_somflow];
% d_3=[eq,zeros(n1,ngap),d_sosvmf,eq-d_sosvmf,zeros(n1,ngap),d_sosvmflow,eq-d_sosvmflow];

%% figure;
inds1=1:500;

dn1=eq-d_mf;
dn2=eq-d_mflow;
dn3=eq-d_somf;
dn4=eq-d_somflow;
dn5=eq-d_sosvmf;
dn6=eq-d_sosvmflow;

figure('units','normalized','Position',[0.0 0.0 .5, 1],'color','w');
% subplot(3,3,1);LO_imagesc(eq(inds1,:),100,2,x,t(inds1)); 
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/2,-0.012,'Raw','color','k','Fontsize',18,'fontweight','bold','HorizontalAlignment','center');
% %text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% %text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/-10,-0.02,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% % text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,1);LO_imagesc(d_mf(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'MF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,2);LO_imagesc(dn1(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.37 0.36],...
    [0.817734959854155 0.756651216011789],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,4,3);LO_imagesc(d_mflow(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'MF+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,4);LO_imagesc(dn2(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.786111111111111 0.776111111111111],...
    [0.829556650246302 0.768472906403935],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,4,5);LO_imagesc(d_somf(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,6);LO_imagesc(dn3(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.371666666666667 0.361666666666667],...
    [0.521185175128967 0.460101431286601],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,4,7);LO_imagesc(d_somflow(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOMF+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(g)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,8);LO_imagesc(dn4(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(h)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.784444444444444 0.774444444444444],...
    [0.526112253209724 0.465028509367358],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,4,9);LO_imagesc(d_sosvmf(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOSVMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(i)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,10);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(j)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.37 0.36],...
    [0.221674876847288 0.160591133004922],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,4,11);LO_imagesc(d_sosvmflow(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOSVMF+LOW (LO)','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(k)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

subplot(3,4,12);LO_imagesc(dn6(inds1,:),100,2,x,t(inds1)); 
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'Removed noise','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.02,'(l)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.783888888888889 0.773888888888889],...
    [0.227588146309293 0.166504402466928],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

rect=[10,10,1];niter=20;eps=0;verb=0;
[simid1]=LO_localsimi(eq-d_mf,d_mf,rect,niter,eps,verb);
[simid2]=LO_localsimi(eq-d_mflow,d_mflow,rect,niter,eps,verb);
[simid3]=LO_localsimi(eq-d_somf,d_somf,rect,niter,eps,verb);
[simid4]=LO_localsimi(eq-d_somflow,d_somflow,rect,niter,eps,verb);
[simid5]=LO_localsimi(eq-d_sosvmf,d_sosvmf,rect,niter,eps,verb);
[simid6]=LO_localsimi(eq-d_sosvmflow,d_sosvmflow,rect,niter,eps,verb);

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(3,2,1);LO_imagesc(simid1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'MF','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.187777777777778 0.174444444444444],...
    [0.846290640394089 0.772413793103448],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

% figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(3,2,2);LO_imagesc(simid2(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'MF+LOW (LO)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.624444444444444 0.611111111111111],...
    [0.843334975369457 0.769458128078816],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,2,3);LO_imagesc(simid3(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOMF','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.188888888888889 0.175555555555556],...
    [0.545798029556649 0.471921182266009],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,2,4);LO_imagesc(simid4(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOMF+LOW (LO)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.627777777777778 0.614444444444444],...
    [0.548753694581278 0.474876847290638],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,2,5);LO_imagesc(simid5(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOSVMF','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.187777777777778 0.174444444444444],...
    [0.245305418719211 0.171428571428571],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);

subplot(3,2,6);LO_imagesc(simid6(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0 1]); colorbar;
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.012,'SOSVMF+LOW (LO)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
%text(n2/2.5,-0.015,'MF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
%text(n2/1.25,-0.015,'MF+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(n2/-10,-0.01,'(f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'arrow',[0.63 0.616666666666667],...
    [0.250231527093595 0.176354679802955],'Color',[1 0 0],'LineWidth',3,...
    'HeadWidth',25,...
    'HeadLength',25);


%%
dn1=eq-d_mf;
dn2=eq-d_mflow;
dn3=eq-d_somf;
dn4=eq-d_somflow;
dn5=eq-d_sosvmf;
dn6=eq-d_sosvmflow;
%% Raw data
t=[0:n1-1]*0.0005;
inds=40;
trace = eq(:,inds);
% length(trace)
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w'); hold on
subplot(7,1,1);
plot(t,trace,'r','linewidth',4); hold on
legend('Raw');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% BP+SOMF+FK method
t=[0:n1-1]*0.0005;
% inds=20;
trace1 = d_mf(:,inds);
%figure;
subplot(7,1,2);
plot(t,trace1,'r','linewidth',4); hold on
legend('MF');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% BP+SOMF+FK method 
% inds=20;
trace2 = d_mflow(:,inds);
%figure;
subplot(7,1,3);
plot(t,trace2,'r','linewidth',4); hold on
legend('MF+LOW (LO)');
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% Proposed LO method
% inds=30;
trace3 = d_somf(:,inds);
subplot(7,1,4)
plot(t,trace3,'r','linewidth',4); hold on
legend('SOMF');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%%
% inds=30;
trace4 = d_somflow(:,inds);
subplot(7,1,5)
plot(t,trace4,'r','linewidth',4); hold on
legend('SOMF+LOW (LO)');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% Proposed LO method
% inds=30;

trace5 = d_sosvmf(:,inds);
subplot(7,1,6)
plot(t,trace5,'r','linewidth',4); hold on
legend('SOSVMF');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(f)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% inds=30;
trace6 = d_sosvmflow(:,inds);
subplot(7,1,7)
plot(t,trace6,'r','linewidth',4); hold on
legend('SOSVMF+LOW (LO)');
% xlim([1 2500])
xlabel('Time (s)')
ylabel('Amplitude')
 ylim([-450 600])
 text(-0.05,1000,'(g)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
%% S/N estimation
a=20*log(rms(trace1))/log(rms(trace))
b=20*log(rms(trace2))/log(rms(trace))
c=20*log(rms(trace3))/log(rms(trace))
d=20*log(rms(trace4))/log(rms(trace))
e=20*log(rms(trace5))/log(rms(trace))
f=20*log(rms(trace6))/log(rms(trace))

% a =
% 
%    17.0100
% 
% 
% b =
% 
%    18.2383
% 
% 
% c =
% 
%    18.0814
% 
% 
% d =
% 
%    18.6213
% 
% 
% e =
% 
%    18.5039
% 
% 
% f =
% 
%    18.8154
