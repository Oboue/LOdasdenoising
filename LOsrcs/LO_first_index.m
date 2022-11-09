function [ i0 ] = amf_first_index( i, j, dim, n, s )
%% Find first index for multidimensional transforms 
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
% i: dimension [0...dim-1]; j: line coordinate; dim: number of dimensions; n: box size [dim]; s: step [dim]
n123 = 1;
i0 = 0;
for k = 0 : dim-1
    if (k == i)
        continue;
    end
	ii = floor(mod((j/n123), n(k+1)));
	n123 = n123 * n(k+1);	
	i0 = i0 + ii * s(k+1);
end

end
