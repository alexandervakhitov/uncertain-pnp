function [M, viz_sigmas] = compute_M_uncertain(U, Alph, A, Sigmas3D, Sigmas2D, Cw, Rest, test, ...
    Alph_S, Alph_E, l2d, Sigmas3DLines, Sigmas2DLines)

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

np = size(Alph,1); %number of 3d points
nl = size(Alph_S,1); %number of 3d points

fu=A(1,1);
fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M = 2 * (np + nl);
ncols_M = 12;
M = zeros(nrows_M,ncols_M);

if (size(Rest, 1) < 3)       
    avg_depth = test;
    Sigmas2D_2D = Sigmas2D * avg_depth * avg_depth;
    Sigmas2D_2D_S = Sigmas2DLines * avg_depth * avg_depth;
    Sigmas2D_2D_E = Sigmas2DLines * avg_depth * avg_depth;
    Sigmas3DCam =  Sigmas3D;
    
    Sigmas3DCamLines = Sigmas3DLines;
else
    Cp_est = Cw*Rest'+repmat(test', 4, 1);
    pts_est = Alph * Cp_est;
    depths_est = pts_est(:, 3);   
    depth_est_rep = repmat(reshape(depths_est, 1, 1, np), 2, 2, 1);
    Sigmas2D_2D = Sigmas2D .* depth_est_rep .* depth_est_rep;
    
    pts_s_est = Alph_S * Cp_est;
    depths_s_est = pts_s_est(:, 3);    
    Sigmas2D_2D_S = Sigmas2DLines .* depths_s_est .* depths_s_est;
    
    pts_e_est = Alph_E * Cp_est;
    depths_e_est = pts_e_est(:, 3);    
    Sigmas2D_2D_E = Sigmas2DLines .* depths_e_est .* depths_e_est;
        
    Sigmas3DCam = zeros(size(Sigmas3D));    
    for i = 1:np
        Sigmas3DCam(:, :, i) = Rest * Sigmas3D(:, :, i) * Rest';
    end
    
    Sigmas3DCamLines = zeros(size(Sigmas3DLines));
    for i = 1:nl
        Sigmas3DCamLines(1:3, 1:3, i) = Rest * Sigmas3DLines(1:3, 1:3, i) * Rest';
        Sigmas3DCamLines(4:6, 4:6, i) = Rest * Sigmas3DLines(4:6, 4:6, i) * Rest';
    end
end

viz_sigmas = zeros(np, 1);

for i = 1:np
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    a4=Alph(i,4);
    
    ui=U(i,1);
    vi=U(i,2);
    
    %generate submatrix M
    M_=[a1*fu, 0, a1*(u0-ui), a2*fu, 0, a2*(u0-ui), a3*fu, 0, a3*(u0-ui), a4*fu, 0, a4*(u0-ui);
        0, a1*fv, a1*(v0-vi), 0, a2*fv, a2*(v0-vi), 0, a3*fv, a3*(v0-vi), 0, a4*fv, a4*(v0-vi)];
    Sigma2D_2D = Sigmas2D_2D(:, :, i);
    Sigma2D_3D = zeros(2,2);
    Sigma3d = Sigmas3DCam(:, :, i);
    du = u0 - ui;
    dv = v0 - vi;
    fu = A(1,1);
    fv = A(2,2);
    Sigma2D_3D(1, 1) = du * du * Sigma3d(3, 3) + 2 * fu * du * Sigma3d(1, 3) + fu * fu * Sigma3d(1, 1);
    Sigma2D_3D(2, 2) = dv * dv * Sigma3d(3, 3) + 2 * fv * dv * Sigma3d(2, 3) + fv * fv * Sigma3d(2, 2);
    Sigma2D_3D(1, 2) = du * dv * Sigma3d(3, 3) + fu * dv * Sigma3d(1, 3) + fv * du * Sigma3d(2, 3) + fu * fv * Sigma3d(1, 2);
    Sigma2D_3D(2, 1) = Sigma2D_3D(1, 2);
    SigmaUnc = Sigma2D_2D+Sigma2D_3D;
    SigmaUncInv = inv(SigmaUnc + 1e-9 * eye(2));    
    SqrtSigmaInv = chol(SigmaUncInv);
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=SqrtSigmaInv*M_;
    
    %visualization
    viz_sigmas(i) = sqrt(trace(SigmaUncInv));
     
end

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
    S6D = Sigmas3DCamLines(:,:,i);    
    Sigma_3d = diag([l*S6D(1:3, 1:3)*l', l*S6D(4:6, 4:6)*l']);
    
    Sigma = Sigma_2d + Sigma_3d + 1e-9*eye(2);
    SigmaInv = inv(Sigma);    
    U = chol(SigmaInv);
    row_ini = 2*np+2*i-1;
    row_end = row_ini+1;
    M(row_ini:row_end,:)=U*M_;
end