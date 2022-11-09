function [ x ] = LO_conjgrad( prec, p, x, dat, niter, ...
                np_conjgrad, nx_conjgrad, nr_conjgrad, nd_conjgrad, ...
                eps_conjgrad, tol_conjgrad, hasp0_conjgrad, r_conjgrad, ...
                sp_conjgrad, gp_conjgrad, sx_conjgrad, gx_conjgrad, ...
                sr_conjgrad, gr_conjgrad, w_weight, n_trianglen, s_trianglen, ...
                nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen,verb )
%% Conjugate-gradient with shaping regularization.
% modified from software Madagascar (University of Texas at Austin)

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
% oper = weight_lop; shape = trianglen_lop; prec is if there is preconditioning

if nargin==26
    verb=1;
end

if (prec)
    d = zeros(nd_conjgrad, 1);
    for i = 0 : nd_conjgrad-1
        d(i+1) = - dat(i+1);
    end
%     prec(false,false,nd_conjgrad,nr_conjgrad,d,r_conjgrad);
else
    for i = 0 : nr_conjgrad-1
        r_conjgrad(i+1) = - dat(i+1);
    end    
end

if (hasp0_conjgrad)
    [p,x] = trianglen_lop(false,false,np_conjgrad,nx_conjgrad,p,x, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen);
    if (prec)
        [x,d] = weight_lop( w_weight, false,false,nx_conjgrad,nd_conjgrad,x,d );
    %     prec(false,true,nd_conjgrad,nr_conjgrad,d,r_conjgrad);    
    else
        [x,r_conjgrad] = weight_lop( w_weight, false,false,nx_conjgrad,nr_conjgrad,x,r_conjgrad );
    end
else
    for i = 0 : np_conjgrad-1
        p(i+1) = 0.0;
    end
    for i = 0 : nx_conjgrad-1
        x(i+1) = 0.0;
    end
end

gnp = 0.0;
dg = 0.0;
g0 = 0.0;

% r0 = sum(r_conjgrad .* r_conjgrad);
r0 = LO_cblas_dsdot(nr_conjgrad,r_conjgrad,1,r_conjgrad,1);


if (r0 == 0.0)
    return
end


% iter = 1

for iter = 0 : niter-1

    for i = 0 : np_conjgrad-1
        gp_conjgrad(i+1) = eps_conjgrad * p(i+1);
    end
    for i = 0 : nx_conjgrad-1
        gx_conjgrad(i+1) = - eps_conjgrad * x(i+1);
    end
    
    if (prec)
        prec(true,true,nd_conjgrad,nr_conjgrad,d,r_conjgrad);  
        [gx_conjgrad,d] = weight_lop( w_weight, true,true,nx_conjgrad,nd_conjgrad,gx_conjgrad,d );
    else
        [gx_conjgrad,r_conjgrad] = LO_weight_lop( w_weight, true,true,nx_conjgrad,nr_conjgrad,gx_conjgrad,r_conjgrad );
    end
    
    [gp_conjgrad,gx_conjgrad] = LO_trianglen_lop(true,true,np_conjgrad,nx_conjgrad,gp_conjgrad,gx_conjgrad, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen);
    [gp_conjgrad,gx_conjgrad] = LO_trianglen_lop(false,false,np_conjgrad,nx_conjgrad,gp_conjgrad,gx_conjgrad, n_trianglen, s_trianglen, nd_trianglen, dim_trianglen, tr_trianglen, tmp_trianglen);
    
    if (prec)
        prec(false,false,nx_conjgrad,nd_conjgrad,gx_conjgrad,d);  
        [d,gr_conjgrad] = LO_weight_lop( w_weight, true,true,nd_conjgrad,nr_conjgrad,d,gr_conjgrad );
    else
        [gx_conjgrad,gr_conjgrad] = LO_weight_lop( w_weight, false,false,nx_conjgrad,nr_conjgrad,gx_conjgrad,gr_conjgrad );
    end    
    
%     gn = sum(gp_conjgrad .* gp_conjgrad);
    gn = LO_cblas_dsdot(np_conjgrad,gp_conjgrad,1,gp_conjgrad,1);
    
    if (iter == 0)
        g0 = gn;
        sp_conjgrad = gp_conjgrad;
        sx_conjgrad = gx_conjgrad;
        sr_conjgrad = gr_conjgrad;
    else
        alpha = gn / gnp;
        dg = gn / g0;
        if ((alpha < tol_conjgrad) || (dg < tol_conjgrad))
            break;
        end
%     sp_conjgrad = gp_conjgrad + alpha * sp_conjgrad;
%     sx_conjgrad = gx_conjgrad + alpha * sx_conjgrad;
%     sr_conjgrad = gr_conjgrad + alpha * sr_conjgrad;
        gp_conjgrad = LO_cblas_saxpy( np_conjgrad, alpha, sp_conjgrad, 1, gp_conjgrad, 1);
        [ gp_conjgrad,sp_conjgrad ] = LO_cblas_sswap( np_conjgrad, gp_conjgrad, 1, sp_conjgrad, 1 );
        
        gx_conjgrad = LO_cblas_saxpy( nx_conjgrad, alpha, sx_conjgrad, 1, gx_conjgrad, 1);
        [ gx_conjgrad,sx_conjgrad ] = LO_cblas_sswap( nx_conjgrad, gx_conjgrad, 1, sx_conjgrad, 1 );
        
        gr_conjgrad = LO_cblas_saxpy( nr_conjgrad, alpha, sr_conjgrad, 1, gr_conjgrad, 1);
        [ gr_conjgrad,sr_conjgrad ] = LO_cblas_sswap( nr_conjgrad, gr_conjgrad, 1, sr_conjgrad, 1 );        

    end
    
    beta = sum(sr_conjgrad .* sr_conjgrad) + eps_conjgrad * (sum(sp_conjgrad .* sp_conjgrad) - sum(sx_conjgrad .* sx_conjgrad));

    if verb
        fprintf('iteration: %d, res: %g !\n',iter,sum(r_conjgrad .* r_conjgrad) / r0);  
    end
%     dg
    
    alpha = - gn / beta;
       
    p = LO_cblas_saxpy( np_conjgrad, alpha, sp_conjgrad, 1, p, 1);
    x = LO_cblas_saxpy( nx_conjgrad, alpha, sx_conjgrad, 1, x, 1); 
    r_conjgrad = LO_cblas_saxpy( nr_conjgrad, alpha, sr_conjgrad, 1, r_conjgrad, 1); 
    
    gnp = gn;
    
    
end




end

