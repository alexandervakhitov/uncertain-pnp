function [M]=compute_M_uncertain_fast(U,Alph,A, sigmas3D, sigmas2D, Cw, Rest, test)

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
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

n=size(Alph,1); %number of 3d points

fu=A(1,1);
fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M=2*n;
ncols_M=12;
M=zeros(nrows_M,ncols_M);
% SM=zeros(nrows_M,ncols_M);


if size(Rest, 1) == 0    
    Sigmas2D_2D = sigmas2D*test*test;%test contains avg point depth
else
    Cp_est = Cw*Rest'+repmat(test', 4, 1);
    est_pts = (Alph * Cp_est)';
    depths = abs(est_pts(3, :))';
    Sigmas2D_2D = sigmas2D.*depths.*depths;
end
alphas = fu*fu*sigmas3D+Sigmas2D_2D;
epsls = sqrt(alphas);
duv = repmat([u0; v0], 1, n) - U';
gamma_norms_sq = sum(duv.^2, 1);
betas = sigmas3D;
delts = epsls ./ gamma_norms_sq .* (-1 + sqrt(1 + betas ./ alphas .* gamma_norms_sq));    
% sigma_est = Sigmas2D_2D + sigmas3D;
% mean_sigma_sqrt = sqrt(1./mean(1./sigma_est));
for i=1:n
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    a4=Alph(i,4);
    
    ui=U(i,1);
    vi=U(i,2);
    
    %generate submatrix M
    M_=[a1*fu, 0, a1*(u0-ui), a2*fu, 0, a2*(u0-ui), a3*fu, 0, a3*(u0-ui), a4*fu, 0, a4*(u0-ui);
        0, a1*fv, a1*(v0-vi), 0, a2*fv, a2*(v0-vi), 0, a3*fv, a3*(v0-vi), 0, a4*fv, a4*(v0-vi)];    
    
%     sigma3d = sigmas3D(i);
%     du = u0 - ui;
%     dv = v0 - vi;
%     fu = A(1,1);
%     fv = A(2,2);
%     
%     duv = [du; dv];        
%     beta = sigmas3D(i);
%     gamma = duv;
%     gamma_norm_sq = gamma' * gamma;
%     delt = epsl / gamma_norm_sq * (-1 + sqrt(1 + beta / alpha * gamma_norm_sq));    

%     SigmaUncSqrt * SigmaUncSqrt - SigmaUnc
    
%     Sigma2D_3D(1, 1) = du * du * sigma3d + fu * fu * sigma3d;
%     Sigma2D_3D(2, 2) = dv * dv * sigma3d + fv * fv * sigma3d;
%     Sigma2D_3D(1, 2) = du * dv * sigma3d;
%     Sigma2D_3D(2, 1) = Sigma2D_3D(1, 2);
%     SigmaUnc = Sigma2D_2D*eye(2)+Sigma2D_3D;
%     SigmaUncSqrt = sqrt(inv(SigmaUnc));%Sigma2D_3D + 

    %fast analytic method    
%     SqrtSigmaInv = 1.0/epsl*(eye(2) - delt / (epsl + delt * gamma_norm_sq) * gamma * gamma');
    
    gamma = duv(:, i);
    SqrtSigmaInv = 1.0/epsls(i)*(eye(2) - delts(i) / (epsls(i) + delts(i) * ...
                   gamma_norms_sq(i)) * gamma * gamma');
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:) = SqrtSigmaInv*M_;    
end

