function [ x,y ] = LO_adjnull( adj,add,nx,ny,x,y )
%% Claerbout-style adjoint zeroing Zeros out the output (unless add is true). 
% Useful first step for any linear operator.
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
% adj : adjoint flag; add: addition flag; nx: size of x; ny: size of y
if(add)
    return
end
if(adj)
    for i = 0 : nx-1
        x(i+1) = 0.0;
    end
else
    for i = 0 : ny-1
        y(i+1) = 0.0;
    end    

end

