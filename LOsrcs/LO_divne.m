function [ rat ] = LO_divne(num, den, rat, eps, niter_divn, n_divn, p_divn, ...
                n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, ...
                tmp_trianglen, np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
                eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, ...
                gp_conjgrad, sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad,verb)
%% N-dimensional smooth division rat=num/den 
% modified from software Madagascar (University of Texas at Austin)
%
% Documentation here !!!!
%%
%   Copyright (C) 2016 Delft University of Technology -- Delphi consortium - Shan Qu
%
%   Modified by Yangkang Chen, Nov 1, 2016
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

if nargin==27
   verb=1;
end

if (eps > 0.0)
    for i = 0 : n_divn-1
        norm = 1.0 / hypot(den(i+1), eps); % sqrt(abs(a).^2 + abs(b).^2)
        num(i+1) = num(i+1) * norm;
        den(i+1) = den(i+1) * norm;
    end
end

% norm = sum(den .* den);
norm = LO_cblas_dsdot(n_divn,den,1,den,1);

if ( norm == 0.0)
    for i = 0 : n_divn-1
        rat(i+1) = 0.0;
    end
    return
end

norm = sqrt(n_divn / norm);

for i = 0 : n_divn-1
    num(i+1) = num(i+1) * norm;
    den(i+1) = den(i+1) * norm;
end

[ w_weight ] = LO_weight_init( den );


[rat ] = LO_conjgrad(false, p_divn, rat, num, niter_divn, np_conjgrad, nx_conjgrad, ...
                nr_conjgrad, nd_conjgrad, eps_conjgrad, tol_conjgrad, hasp0_conjgrad, ...
                r_conjgrad, sp_conjgrad, gp_conjgrad, sx_conjgrad, gx_conjgrad, ...
                sr_conjgrad, gr_conjgrad, w_weight, n_trianglen, s_trianglen, ...
                nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen, verb);

% rat1 .* w_weight - num
end

