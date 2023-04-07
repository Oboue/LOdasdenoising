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
ngap0=1000;

ngap=50;
x=1:n2*3+2*ngap;
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
ns=14;                        % spray radius
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
w=0.03;                       % half width (in percentage) of the cone filter (i.e., w*nk=nwidth)
%
tic
d3=LO_bandpasssosvmffk(din,dt,flo,fhi,nplo,nphi,phase,verb0,niter,liter,order1,eps_dv,eps_cg,tol_cg,rect,verb1,adj,add,n1,n2,ns,order2,eps,ndn,nds,type_mf,ifsmooth,w);
toc
%
%% Denoising using the BP+SOSVMF+FK+curvelet method 
%  Parameter tuning：add the ke              y parameters of the curvelet method

c1=1;                % Type of the transform(0: complex-valued curvelets,1: real-valued curvelets)
c2=1;                % Chooses one of two possibilities for the coefficients at the finest level(1: curvelets,2: wavelets)
c3=0.50;             % Thresholding parameter (alpha)
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
size(din)
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

figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(1,1,1);LO_imagesc(din,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/0.65,-0.025,'Raw DAS data','color','k','Fontsize',40,'fontweight','bold','HorizontalAlignment','center');
% text(-200,-0.05,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.651666666666666 0.693888888888889],...
    [0.743262081559298 0.674432189720557],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-frequency noise'},...
    'LineWidth',5,...
    'HeadWidth',40,...
    'HeadLength',25,...
    'FontWeight','bold',...
    'FontSize',40);
annotation(gcf,'textarrow',[0.733333333333333 0.742777777777778],...
    [0.5366395745477 0.602955665024629],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Signal'},...
    'LineWidth',5,...
    'HeadWidth',40,...
    'HeadLength',25,...
    'FontWeight','bold',...
    'FontSize',40);
annotation(gcf,'textarrow',[0.609444444444444 0.634444444444444],...
    [0.267535043998961 0.349753694581281],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',5,...
    'HeadWidth',40,...
    'HeadLength',25,...
    'FontWeight','bold',...
    'FontSize',40);
annotation(gcf,'textarrow',[0.39111111111111 0.391111111111111],...
    [0.227932235729295 0.178325123152708],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Random','background noise'},...
    'LineWidth',5,...
    'HeadWidth',40,...
    'HeadLength',25,...
    'FontWeight','bold',...
    'FontSize',40);
annotation(gcf,'textarrow',[0.215555555555556 0.195],...
    [0.317920362725042 0.451231527093596],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'   Horizontal noise'},...
    'LineWidth',5,...
    'HeadWidth',40,...
    'HeadLength',25,...
    'FontWeight','bold',...
    'FontSize',40);
%% plot the waveforms
figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,3,1);LO_imagesc(din,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'Raw','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.255555555555555 0.256666666666667],...
    [0.858474107650487 0.835629921259841],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'Random','background noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.243333333333333 0.285555555555556],...
    [0.794493608652902 0.725663716814161],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-frequency noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.296666666666666 0.313333333333333],...
    [0.661623713949702 0.63779527559055],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);
annotation(gcf,'textarrow',[0.25 0.25],...
    [0.609545978488593 0.625984251968504],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'   Horizontal noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(2,3,2);LO_imagesc(d1,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.573333333333333 0.591111111111112],...
    [0.641938674579624 0.615232443125618],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'High-amplitude','erratic noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(2,3,3);LO_imagesc(d2,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP+SOSVMF','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.802222222222222 0.803333333333333],...
    [0.659742828882294 0.632512315270936],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'   Horizontal noise'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(2,3,4);LO_imagesc(d3,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP+SOSVMF+FK','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.187666666666666 0.323698430564284 0.106777777777778 0.0741839762611269],...
    'Color',[1 0 0],...
    'String',{'Random background noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
subplot(2,3,5);LO_imagesc(d4,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.54 0.577777777777778],...
    [0.347227785005334 0.279802955665024],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'     Weaker','    signal energy'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

subplot(2,3,6);LO_imagesc(d5,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.05,'(f)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.812222222222222 0.86],...
    [0.352181179439951 0.277008083495343],'Color',[1 0 0],'TextColor',[1 0 0],...
    'String',{'  Improved','    signal energy'},...
    'LineWidth',2,...
    'HeadWidth',15,...
    'HeadLength',15,...
    'FontWeight','bold',...
    'FontSize',18);

dn1=din-d1;
dn2=d1-d2;
dn3=d2-d3;
dn4=d3-d4;
dn5=dn1+dn2+dn3+dn4;
dn6=din-d5;

figure('units','normalized','Position',[0.0 0.0 0.5 1],'color','w');
subplot(2,3,1);LO_imagesc(dn1,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.025,'(a)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.141111111111111 0.727362204724409 0.206666666666667 0.104155928784759],...
    'Color',[1 0 0],...
    'String',{'','High-frequency','noise'},...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% figure('units','normalized','Position',[0.0 0.0 1 1],'color','w');
subplot(2,3,2);LO_imagesc(dn2,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(b)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.433333333333333 0.747815018812295 0.184444444444445 0.0394088669950738],...
    'Color',[1 0 0],...
    'String','Erratic noise',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off',...
    'EdgeColor','none');

subplot(2,3,3);LO_imagesc(dn3,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(c)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.7 0.751459175361687 0.203333333333334 0.0394088669950712],...
    'Color',[1 0 0],...
    'String','Horizontal noise',...
    'FontWeight','bold',...
    'FontSize',22,...
    'FitBoxToText','off',...
    'EdgeColor','none');

subplot(2,3,4);LO_imagesc(dn4,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(d)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.138888888888889 0.210115976882199 0.202222222222222 0.09949261083742],...
    'Color',[1 0 0],...
    'String',{'Random','background noise'},...
    'FontWeight','bold',...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');

subplot(2,3,5);LO_imagesc(dn5,100,2,x,t);caxis([-25,25]); hold on
% subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/0.65,-0.025,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
text(-500,-0.05,'(e)','color','k','Fontsize',25,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textbox',...
    [0.41 0.40255905511811 0.234444444444445 0.0539590783910575],...
    'Color',[1 0 0],...
    'String',{'','High-frequency noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.495555555555555 0.363188976377952 0.163333333333333 0.0826771653543306],...
    'Color',[1 0 0],...
    'String',{'','+'},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.434444444444444 0.309055118110236 0.184444444444445 0.0756097125790205],...
    'Color',[1 0 0],...
    'String',{'','High-amplitude','   Erratic noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.494444444444444 0.27263779527559 0.163333333333333 0.0826771653543305],...
    'Color',[1 0 0],...
    'String',{'','+'},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.424444444444445 0.237679647802632 0.173333333333333 0.0394088669950713],...
    'Color',[1 0 0],...
    'String','Horizontal noise',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.495555555555555 0.194881889763778 0.161111111111112 0.0826771653543306],...
    'Color',[1 0 0],...
    'String',{'','+'},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(gcf,'textbox',...
    [0.421111111111111 0.138779527559055 0.202222222222222 0.0527188239400912],...
    'Color',[1 0 0],...
    'String',{'       Random','background noise'},...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'EdgeColor','none');
% %%
% inds1=800:1000;
% 
% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(2,3,1);LO_imagesc(din(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',11,'fontweight','bold');
% xlabel('Channel','Fontsize',11,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',11,'Fontweight','bold');
% text(n2/0.6,0.397,'Raw','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% annotation(gcf,'textbox',...
%     [0.182222222222222 0.642574257425743 0.117777777777778 0.0762376237623769],...
%     'Color',[1 0 0],...
%     'String',{'Mixture od strong noise'},...
%     'FontWeight','bold',...
%     'FontSize',20,...
%     'FitBoxToText','off');
% 
% subplot(2,3,2);LO_imagesc(d1(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',11,'fontweight','bold');
% xlabel('Channel','Fontsize',11,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',11,'Fontweight','bold');
% annotation(gcf,'textarrow',[0.567777777777777 0.591111111111111],...
%     [0.660945386197303 0.628712871287129],'Color',[1 0 0],'TextColor',[1 0 0],...
%     'String',{'High-amplitude','erratic noise'},...
%     'LineWidth',2,...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FontName','Helvetica Neue');
% text(n2/0.6,0.397,'BP','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,3);LO_imagesc(d2(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% annotation(gcf,'textarrow',[0.807222222222221 0.797222222222222],...
%     [0.718380941447575 0.756371488288571],'Color',[1 0 0],'TextColor',[1 0 0],...
%     'String',{'Horizontal noise'},...
%     'LineWidth',2,...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FontName','Helvetica Neue');
% text(n2/0.6,0.397,'BP+SOSVMF','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,4);LO_imagesc(d3(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% annotation(gcf,'textbox',...
%     [0.171555555555555 0.16453328950338 0.124 0.0810000000000005],...
%     'Color',[1 0 0],...
%     'String',{'Random','background','noise'},...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FitBoxToText','off',...
%     'EdgeColor','none');
% text(n2/0.6,0.397,'BP+SOSVMF+FK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,5);LO_imagesc(d4(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% annotation(gcf,'textarrow',[0.56 0.57],...
%     [0.210817232648548 0.170807779489544],'Color',[1 0 0],'TextColor',[1 0 0],...
%     'String',{'Weaker signal','energy'},...
%     'LineWidth',2,...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FontName','Helvetica Neue');
% text(n2/0.6,0.397,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,3,6);LO_imagesc(d5(inds1,:),100,2,x,t(inds1));caxis([-25,25]); hold on
% ylabel('Time (s)','Fontsize',11,'fontweight','bold');
% xlabel('Channel','Fontsize',11,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',11,'Fontweight','bold');
% annotation(gcf,'textarrow',[0.833333333333333 0.848333333333333],...
%     [0.211839978431551 0.173439048562934],'Color',[1 0 0],'TextColor',[1 0 0],...
%     'String',{'Stronger signal','energy'},...
%     'LineWidth',2,...
%     'HorizontalAlignment','center',...
%     'FontWeight','bold',...
%     'FontSize',18,...
%     'FontName','Helvetica Neue');
% text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');

% text(-200,0.25,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(n2/3.5,0.255,'Raw data','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.65,0.255,'BP','color','k','Fontsize',12,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.38,0.255,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-200,0.395,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% text(n2/2.0,0.395,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.65,0.395,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(n2/0.38,0.395,'LO','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% print(gcf,'-depsc','-r300','fig2.eps');
% %%
% % h
% rect=[30,30,1];niter=20;eps=0;verb=0;
% [simi1]=LO_localsimi(din-d1,d1,rect,niter,eps,verb);
% [simi2]=LO_localsimi(din-d2,d2,rect,niter,eps,verb);
% [simi3]=LO_localsimi(din-d3,d3,rect,niter,eps,verb);
% [simi4]=LO_localsimi(din-d4,d4,rect,niter,eps,verb);
% [simi5]=LO_localsimi(din-d5,d5,rect,niter,eps,verb);
% 
% d_sim1=[simi1,zeros(n1,ngap),simi2,zeros(n1,ngap),simi3];
% d_sim2=[simi3,zeros(n1,ngap),simi5];
% 
% inds1=800:1000;
% 
% dn1=din-d1;
% dn2=din-d2;
% dn3=din-d3;
% dn4=din-d4;
% dn5=din-d5;
% 
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,5,1);LO_imagesc(dn1(inds1,:),100,2,x,t(inds1));caxis([-25,25]); 
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(a)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,2);LO_imagesc(dn2(inds1,:),100,2,x,t(inds1));caxis([-25,25]);
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(b)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,3);LO_imagesc(dn3(inds1,:),100,2,x,t(inds1));caxis([-25,25]);
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF+FK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(c)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,4);LO_imagesc(dn4(inds1,:),100,2,x,t(inds1));caxis([-25,25]);
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(d)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,5);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));caxis([-25,25]); 
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(e)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,6,6);LO_imagesc(dn5(inds1,:),100,2,x,t(inds1));caxis([-25,25]); 
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% %
% figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(2,5,6);LO_imagesc(simi1(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(f)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,7);LO_imagesc(simi2(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(g)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,8);LO_imagesc(simi3(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF+FK','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(h)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(2,5,9);LO_imagesc(simi4(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.6,0.397,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(i)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% subplot(2,5,10);LO_imagesc(simi5(inds1,:),100,1,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% text(n2/0.65,0.397,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',11,'fontweight','bold','HorizontalAlignment','center');
% text(-300,0.393,'(j)','color','k','Fontsize',16,'fontweight','bold','HorizontalAlignment','center');
% 
% 
% % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
% subplot(3,2,1);LO_imagesc(simi1(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',10,'fontweight','bold');
% xlabel('Channel','Fontsize',10,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
% text(n2/0.7,-0.04,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% caxis([0,0.25]);
% 
% subplot(3,2,2);LO_imagesc(simi2(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% xlabel('Channel','Fontsize',14,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% text(n2/0.7,-0.04,'BP+SOSVMF','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% caxis([0,0.25]);
% 
% subplot(3,2,3);LO_imagesc(simi3(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% xlabel('Channel','Fontsize',14,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% text(n2/0.7,-0.04,'BP+SOSVMF+FK','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% caxis([0,0.25]);
% 
% % figure('units','normalized','Position',[0.0 0.0 1, 1],'color','w');
% subplot(3,2,4);LO_imagesc(simi4(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% xlabel('Channel','Fontsize',14,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% text(n2/0.7,-0.04,'BP+SOSVMF+FK+Curvelet','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% caxis([0,0.25]);
% 
% subplot(3,2,5);LO_imagesc(simi5(inds1,:),95,1,x,t);colormap(jet);caxis([0,0.25]); hold on
% % subplot(2,3,1);LO_imagesc(simi1(inds1,:),100,2,x,t(inds1));colormap(jet);caxis([0,0.25]); hold on
% ylabel('Time (s)','Fontsize',14,'fontweight','bold');
% xlabel('Channel','Fontsize',14,'fontweight','bold');
% set(gca,'Linewidth',2,'Fontsize',14,'Fontweight','bold');
% text(n2/0.65,-0.04,'BP+SOSVMF+FK+Curvelet+LOW (LO)','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');
% text(-300,-0.1,'(e)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
% c = colorbar;c.Label.String = 'Local similarity';c.Label.FontSize = 10;%c.Label.FontWeight = bold;
% caxis([0,0.25]);
% 
% print(gcf,'-depsc','-r200','fig10.eps');
