function [ x, triangle_struct ] = LO_smooth2( triangle_struct, o, d, x, der)
%% apply triangle smoothing 
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
% triangle_struct: smoothing object; o and d:trace sampling ; der: if derivative; x: data
triangle_struct.tmp = LO_triple2(o, d, triangle_struct.nx, triangle_struct.nb, x, triangle_struct.tmp, triangle_struct.box, triangle_struct.wt);
triangle_struct.tmp = LO_doubint2(triangle_struct.np, triangle_struct.tmp, (triangle_struct.box || der));
x = LO_fold2(o, d, triangle_struct.nx, triangle_struct.nb, triangle_struct.np, x, triangle_struct.tmp);

end

