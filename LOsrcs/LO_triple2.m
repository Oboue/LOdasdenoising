function [ tmp ] = LO_triple2( o, d, nx, nb, x, tmp, box, wt )
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
for i = 0 : nx+2*nb-1
    tmp(i+1) = 0;
end
if (box)
    tmp(1+1:end)     = LO_cblas_saxpy(nx,  +wt,x(1+o:end),d,tmp(1+1:end),   1); % y += a*x
	tmp(1+2*nb:end)  = LO_cblas_saxpy(nx,  -wt,x(1+o:end),d,tmp(1+2*nb:end),1);
else
    tmp              = LO_cblas_saxpy(nx,  -wt,x(1+o:end),d,tmp,            1); % y += a*x
	tmp(1+nb:end)    = LO_cblas_saxpy(nx,2.*wt,x(1+o:end),d,tmp(1+nb:end),  1);
	tmp(1+2*nb:end)  = LO_cblas_saxpy(nx,  -wt,x(1+o:end),d,tmp(1+2*nb:end),1);
end



end
