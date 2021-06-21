function M=compute_M_uncertain_alph(U,Alph,A, Sigmas3DAlph, Sigmas2D, Cw)

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


Cp_est = Cw*Rest'+repmat(test', 4, 1);

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
    
    est_pt = ([a1 a2 a3 a4] * Cp_est)';
    
    Sigma2D_2D = Sigmas2D(:,:,i)*est_pt(3)*est_pt(3);
%     Sigma2d(0, 0) = du * du * Sigma3d(2, 2) + 2 * fu * du * Sigma3d(0, 2) + fu * fu * Sigma3d(0, 0);
%     Sigma2d(1, 1) = dv * dv * Sigma3d(2, 2) + 2 * fv * dv * Sigma3d(1, 2) + fv * fv * Sigma3d(1, 1);
%     Sigma2d(0, 1) = du * dv * Sigma3d(2, 2) + fu * dv * Sigma3d(0, 2) + fv * du * Sigma3d(1, 2) +
%                         fu * fv * Sigma3d(0, 1);
%     Sigma2d(1, 0) = Sigma2d(0, 1);
    Sigma2D_3D = zeros(2,2);
    Sigma3d = Rest * Sigmas3D(:,:,i)*Rest';
    du = u0 - ui;
    dv = v0 - vi;
    fu = A(1,1);
    fv = A(2,2);
    Sigma2D_3D(1, 1) = du * du * Sigma3d(3, 3) + 2 * fu * du * Sigma3d(1, 3) + fu * fu * Sigma3d(1, 1);
    Sigma2D_3D(2, 2) = dv * dv * Sigma3d(3, 3) + 2 * fv * dv * Sigma3d(2, 3) + fv * fv * Sigma3d(2, 2);
    Sigma2D_3D(1, 2) = du * dv * Sigma3d(3, 3) + fu * dv * Sigma3d(1, 3) + fv * du * Sigma3d(2, 3) + fu * fv * Sigma3d(1, 2);
    Sigma2D_3D(2, 1) = Sigma2D_3D(1, 2);
    SigmaUnc = Sigma2D_2D+Sigma2D_3D;
    SigmaInv = inv(SigmaUnc);%Sigma2D_3D + 
    
    [Us,Ss,Vs] = svd(SigmaInv);
    SqrtSigmaInv = sqrt(Ss)*Vs';
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=SqrtSigmaInv*M_;
     
end

