function [ np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, eps_conjgrad, ....
           tol_conjgrad, hasp0_conjgrad, r_conjgrad, sp_conjgrad, gp_conjgrad, ...
           sx_conjgrad, gx_conjgrad, sr_conjgrad, gr_conjgrad ] = ...
           LO_conjgrad_init( np1, nx1, nd1, nr1, eps1, tol1, hasp01 )
%% initialize Conjugate-gradient with shaping regularization
% modified from software Madagascar (University of Texas at Austin)

%%
%   Copyright (C) 2016 --   Shan Qu
%                           Delphi consortium 
%                           Delft University of Technology 
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
% np1: preconditioned size; nx1: model size; nd1: data size; nr1: residual size; eps1: scaling;
% tol1: tolerance; hasp01: if has initial model 

np_conjgrad = np1; 
nx_conjgrad = nx1;
nr_conjgrad = nr1;
nd_conjgrad = nd1;
eps_conjgrad = eps1*eps1;
tol_conjgrad = tol1;
hasp0_conjgrad = hasp01;

r_conjgrad = zeros(nr_conjgrad, 1);  
sp_conjgrad = zeros(np_conjgrad, 1);
gp_conjgrad = zeros(np_conjgrad, 1);
sx_conjgrad = zeros(nx_conjgrad, 1);
gx_conjgrad = zeros(nx_conjgrad, 1);
sr_conjgrad = zeros(nr_conjgrad, 1);
gr_conjgrad = zeros(nr_conjgrad, 1);

end

