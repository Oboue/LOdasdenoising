function [ x ] = LO_fold2(o, d, nx, nb, np, x, tmp)
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
% copy middle 
for i = 0 : nx-1
    x(o+i*d+1) = tmp(i+nb+1);
end

% reflections from the right side 
for j = nb+nx : nx : np
    if (nx <= np-j)
        for i = 0 : nx-1
            x(o+(nx-1-i)*d+1) = x(o+(nx-1-i)*d+1) + tmp(j+i+1);
        end
    else
        for i = 0 : np-j-1
            x(o+(nx-1-i)*d+1) = x(o+(nx-1-i)*d+1) + tmp(j+i+1);
        end
    end
    j = j + nx;
    if (nx <= np-j)
        for i = 0 : nx-1
            x(o+i*d+1) = x(o+i*d+1) + tmp(j+i+1);
        end

    else
        for i = 0 : np-j-1
            x(o+i*d+1) = x(o+i*d+1) + tmp(j+i+1);
        end
    end
end
% 
%     reflections from the left side
for j = nb : -nx : 0
    if (nx <= j)
        for i = 0 : nx-1
            x(o+i*d+1) = x(o+i*d+1) + tmp(j-1-i+1);
        end
    else
        for i = 0 : j-1
            x(o+i*d+1) = x(o+i*d+1) + tmp(j-1-i+1);
        end
    end
    j = j - nx;
    if (nx <= j)
        for i = 0 : nx-1
            x(o+(nx-1-i)*d+1) = x(o+(nx-1-i)*d+1) + tmp(j-1-i+1);
        end
    else
        for i = 0 : j-1
            x(o+(nx-1-i)*d+1) = x(o+(nx-1-i)*d+1) + tmp(j-1-i+1);
        end
    end
end


end

