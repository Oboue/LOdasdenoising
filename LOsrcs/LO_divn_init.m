function [ niter_divn, n_divn, p_divn, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad] = LO_divn_init( ndim, nd, ndat, nbox, niter1 )
%% initialize N-dimensional smooth division
% modified from software Madagascar (University of Texas at Austin)
%%
%   Copyright (C) 2016 Delft University of Technology -- Delphi consortium - Shan Qu
%   
%   Slightly modified by Yangkang Chen, Nov 1, 2016
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%   
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%%

niter_divn = niter1;
n_divn = nd;

[n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen ] = LO_trianglen_init( ndim, nbox, ndat );

[ np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, eps_conjgrad, tol_conjgrad, ...
    hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, sx_conjgrad, gx_conjgrad, ...
    sr_conjgrad, gr_conjgrad ] = LO_conjgrad_init(nd, nd, nd, nd, 1.0, 1e-6, false);
p_divn = zeros(nd, 1);

end

