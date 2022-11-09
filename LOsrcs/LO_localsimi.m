function [ simi ] = LO_localsimi(d1,d2,rect,niter,eps, verb)
%  LO_LOCALSIMI: calculate local similarity between two datasets
%
%  IN   d1:   	input data 1
%       d2:     input data 2
%       verb:   verbosity flag (default: 0)
%
%  OUT  simi:  	calculated local similarity, which is of the same size as d1 and d2
%
%  Copyright (C) 2016 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  Reference:   
%  				1. Chen, Y. and S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, , 80, WD1-WD9. (Note that when local similarity is used for noise suppression purposes, this reference must be cited.)
%               2. Local seismic attributes, Fomel, Geophysics, 2007
%
% DEMO
% test/test_localsimi.m 

addpath(genpath('~/chenyk.data/dip_estimation_fomel_matlab'));

if nargin<2
    error('Input data 1 and data 2 must be provided!');
end

[n1,n2,n3]=size(d1);

if nargin==2
    rect=ones(3,1);
    if(n1==1) error('data must be a vector or a matrix!');
    else
        rect(1)=20;
    end
    if(n2~=1) rect(2)=10;end
    if(n3~=1) rect(3)=10;end
    niter=50;
    eps=0.0;
    verb=1;
end;

if nargin==3
   niter=50;
   eps=0.0;
   verb=1;
end

if nargin==4
   eps=0.0;
   verb=1;
end

if nargin==5
   verb=1;
end

%eps=0.0;

nd=n1*n2*n3;
n_dat=[n1,n2,n3];

[ niter_divn, n_divn, p_divn, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, ...
    tr_trianglen, tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
    eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, ...
    sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad] ...
    = LO_divn_init(3, nd, n_dat, rect, niter);

ratio=zeros(size(d1));

[ ratio ] = LO_divne( d2, d1, ratio, eps, ...
    niter_divn, n_divn, p_divn, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, ...
    tr_trianglen, tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
    eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, ...
    sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad, verb);

%%% opposite direction
[ niter_divn, n_divn, p_divn, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, ...
    tr_trianglen, tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
    eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, ...
    sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad] ...
    = LO_divn_init(3, nd, n_dat, rect, niter);

ratio1=zeros(size(d1));
[ ratio1 ] = LO_divne( d1, d2, ratio1, eps, ...
    niter_divn, n_divn, p_divn, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, ...
    tr_trianglen, tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
    eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, ...
    sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad, verb);

simi=sqrt(abs(ratio.*ratio1));






