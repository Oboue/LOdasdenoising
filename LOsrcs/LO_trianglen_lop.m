function [ x, y ] = LO_trianglen_lop( adj, add, nx, ny, x, y, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen)
%% N-D triangle smoothing as a linear operator
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
[ x,y ] = LO_adjnull (adj,add,nx,ny,x,y);

if (adj)
    for i = 0 : nd_trianglen-1
        tmp_trianglen(i+1) = y(i+1);
    end
else
    for i = 0 : nd_trianglen-1
        tmp_trianglen(i+1) = x(i+1);
    end  
end

for i = 0 : dim_trianglen-1
    if not(isempty(tr_trianglen{i+1}))
        for j = 0 : (nd_trianglen / n_trianglen(i+1) - 1)
            i0 = LO_first_index(i,j,dim_trianglen,n_trianglen,s_trianglen);
            [tmp_trianglen, tr_trianglen{i+1}] = LO_smooth2( tr_trianglen{i+1}, i0, s_trianglen(i+1), tmp_trianglen,false);
        end
    end
end

if (adj)
    for i = 0 : nd_trianglen-1
        x(i+1) = x(i+1) + tmp_trianglen(i+1);
    end
else
    for i = 0 : nd_trianglen-1
        y(i+1) = y(i+1) + tmp_trianglen(i+1);
    end
end



end

