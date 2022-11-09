function [n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen ] = LO_trianglen_init( ndim, nbox, ndat )
%% initialize N-D triangle smoothing as a linear operator
% modified from software Madagascar (University of Texas at Austin)

%%
%   Copyright (C) 2016 Delft University of Technology -- Delphi consortium - Shan Qu
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
% nbox: triangle radius [ndim]; ndat: data dimensions [ndim] 

dim_trianglen = ndim;
n_trianglen = zeros(dim_trianglen, 1);
tr_trianglen = cell(dim_trianglen,1);

nd_trianglen = 1;
s_trianglen = zeros(dim_trianglen, 1);

for i = 0 : dim_trianglen-1
    if (nbox(i+1) > 1)
        tr_trianglen{i+1} = LO_triangle_init( nbox(i+1), ndat(i+1), false);
    end
    s_trianglen(i+1) = nd_trianglen;
    n_trianglen(i+1) = ndat(i+1);
    nd_trianglen = nd_trianglen * ndat(i+1);
end

tmp_trianglen = zeros(nd_trianglen, 1);
        

end
