function [Xc_opt,R_opt,T_opt,err_opt,iter,Cc,errs]=optimize_betas_gauss_newton(Km,...
    Cw,Beta0,Xs,Xe,xs,xe,A,Xw, U,Alpha, Alpha_S, ptMean, RA, ma)

% COMPUTE_BETAS_GAUSS_NEWTON  
%
%       Km: vector of the kernel
%       Cw: position of the control point in world coordinates
%       Beta0: initial guess of the betas
%
%
% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, October 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 


n=size(Beta0,2);

[Beta_opt,err,iter] = epnp_lines.gauss_newton(Km,Cw,Beta0);

%Extract control point camera coordinates from Betas and Kernel
X=zeros(12,1);
for i=1:n
   X=X+Beta_opt(i)*Km(:,i); 
end

Cc=zeros(4,3);
for i=1:4
    Cc(i,:)=X(3*i-2:3*i);
end


%check sign of the determinant (keep orienation of the control points)
s_Cw = epnp_lines.sign_determinant(Cw);
s_Cc = epnp_lines.sign_determinant(Cc);
Cc=Cc*(s_Cw/s_Cc);

%Reconstruct and compute error=

if (length(Xw) > 3)
    Xc_opt = Cc'*Alpha; %reconstruction: points in camera coordinate system
    [R_opt,T_opt] = epnp_lines.getrotT(ptMean, RA, Cc, ma);  %solve exterior orientation
    [err_opt,Urep_opt] = epnp_lines.reprojection_error_usingRT(Xw,U,R_opt,T_opt,A);
    errs = [];
else    
    Xc_opt = Cc'*Alpha_S;
    [R_opt,T_opt] = epnp_lines.getrotT(ptMean, RA, Cc, ma); %solve exterior orientation
    [err_opt,~,~,errs] = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R_opt,T_opt,A);
    errs = [];
end

end




