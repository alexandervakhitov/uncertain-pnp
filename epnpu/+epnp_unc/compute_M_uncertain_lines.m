function M=compute_M_uncertain_lines(U,Alph,A, sigmas3D, sigmas2D, Cw, Rest, ...
    test, Alph_S, Alph_E, l2d, Sigmas3DLines, Sigmas2DLines)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       A. Vakhitov, L. Ferraz, A. Agudo, F. Moreno-Noguer
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       CVPR 2021
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
% Alexander Vakhitov, alexander.vakhitov@gmail.com, SLAMcore, Jun. 2021
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 
 

np=size(Alph,1); %number of 3d points
nl = size(Alph_S, 1);

fu=A(1,1);
fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M=2*(np+nl);
ncols_M=12;
M=zeros(nrows_M,ncols_M);

if size(Rest, 1) == 0    
    Sigmas2D_2D = sigmas2D*test*test;%test contains avg point depth
    Sigmas2D_2D_S = Sigmas2DLines*test*test;
    Sigmas2D_2D_E = Sigmas2DLines*test*test;
else
    Cp_est = Cw*Rest'+repmat(test', 4, 1);
    est_pts = (Alph * Cp_est)';
    depths = abs(est_pts(3, :))';
    Sigmas2D_2D = sigmas2D.*depths.*depths;
    
    start_est_pts = (Alph_S * Cp_est)';
    end_est_pts = (Alph_E * Cp_est)';
    start_depths = abs(start_est_pts(3, :))';
    end_depths = abs(end_est_pts(3, :))';
    Sigmas2D_2D_S = Sigmas2DLines.*start_depths.*start_depths;
    Sigmas2D_2D_E = Sigmas2DLines.*end_depths.*end_depths;
end
alphas = fu*fu*sigmas3D+Sigmas2D_2D + 1e-14;
epsls = sqrt(alphas);
duv = repmat([u0; v0], 1, np) - U';
gamma_norms_sq = sum(duv.^2, 1);
betas = sigmas3D;
delts = epsls ./ gamma_norms_sq .* (-1 + sqrt(1 + betas ./ alphas .* gamma_norms_sq));    
for i=1:np
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    a4=Alph(i,4);
    
    ui=U(i, 1);
    vi=U(i, 2);
    
    %generate submatrix M
    M_=[a1*fu, 0, a1*(u0-ui), a2*fu, 0, a2*(u0-ui), a3*fu, 0, a3*(u0-ui), a4*fu, 0, a4*(u0-ui);
        0, a1*fv, a1*(v0-vi), 0, a2*fv, a2*(v0-vi), 0, a3*fv, a3*(v0-vi), 0, a4*fv, a4*(v0-vi)];    
    
    gamma = duv(:, i);
    SqrtSigmaInv = 1.0/epsls(i)*(eye(2) - delts(i) / (epsls(i) + delts(i) * ...
                   gamma_norms_sq(i)) * gamma * gamma');
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:) = SqrtSigmaInv*M_;    
end

l_norms_sq = sum(l2d.^2, 1);

for i = 1:nl
    l = l2d(:, i)';
    a1s=Alph_S(i,1);
    a2s=Alph_S(i,2);
    a3s=Alph_S(i,3);
    a4s=Alph_S(i,4);

    a1e=Alph_E(i,1);
    a2e=Alph_E(i,2);
    a3e=Alph_E(i,3);
    a4e=Alph_E(i,4);
    M_ = [ a1s*l, a2s*l, a3s*l, a4s*l;
           a1e*l, a2e*l, a3e*l, a4e*l];

    Sigma_2d = diag([Sigmas2D_2D_S(i), Sigmas2D_2D_E(i)]);
    
    Sigma_3d = eye(2) * l_norms_sq(i) * Sigmas3DLines(i);
    Sigma = Sigma_2d + Sigma_3d + 1e-9*eye(2);
    SigmaInv = inv(Sigma);
    
    U = chol(SigmaInv);
    row_ini = 2*np+2*i-1;
    row_end = row_ini+1;
    M(row_ini:row_end,:)=U*M_;
end
