function Km = compute_M_ver3_unc(U, Alph,A, Alph_S, Alph_E, l2d, Sigmas3D, Sigmas2D, SigmasLines, SigmasLines2D, Cw, Rest, test)

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

n=size(Alph, 2); %number of 3d points
nl = size(Alph_S, 2);

fu = A(1,1);
fv = A(2,2);
u0 = A(1,3);
v0 = A(2,3);

nrows_M = 2 * (n+nl);
nrows_M = 2 * (n);
ncols_M = 12;

At = Alph';
Ut = U';
M1 = [fu*ones(n, 4) zeros(n, 4) repmat(u0-Ut(:, 1), 1, 4)].*repmat(At, 1, 3);
M2 = [zeros(n, 4) fv*ones(n, 4) repmat(v0-Ut(:, 2), 1, 4)].*repmat(At, 1, 3);

Mpts = zeros(2*size(M1,1), size(M1,2));
Mpts(1:2:end,:) = M1;
Mpts(2:2:end,:) = M2;

Cp_est = Cw*Rest'+repmat(test', 4, 1);
Alpht = Alph';

for i = 1:n
    a1=Alpht(i,1);
    a2=Alpht(i,2);
    a3=Alpht(i,3);
    a4=Alpht(i,4);
    est_pt = ([a1 a2 a3 a4] * Cp_est)';
    
    ui=U(1,i);
    vi=U(2,i);

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
    
    Mpts(row_ini:row_end, :) = SqrtSigmaInv*Mpts(row_ini:row_end, :);
end

%Me1 = M1'*M1 + M2'*M2;
Me1 = Mpts'*Mpts;
%from here works only for u0 = v0 = 0, fu = fv

AS = [Alph_S Alph_E]';
l22 = repmat(l2d', 2, 1);
Me = [repmat(l22(:,1),1,4) repmat(l22(:,2),1,4) repmat(l22(:,3),1,4) ] .* repmat(AS, 1, 3);

Mlines = zeros(2*nl, 12);
Mlines(1:2:end, :) = Me(1:nl,:);
Mlines(2:2:end, :) = Me(nl+1:end,:);

Xsc = (Alph_S' * Cp_est)';
Xec = (Alph_E' * Cp_est)';

for i = 1:nl
    Sigma2 = zeros(2,2);
    Sigma2(1,1) = l2d(:, i)'*Rest*SigmasLines(1:3,1:3,i)*Rest'*l2d(:, i);
    Sigma2(1,2) = l2d(:, i)'*Rest*SigmasLines(1:3,4:6,i)*Rest'*l2d(:, i);
    Sigma2(2,2) = l2d(:, i)'*Rest*SigmasLines(4:6,4:6,i)*Rest'*l2d(:, i);        
    Sigma2(2,1) = l2d(:, i)'*Rest*SigmasLines(4:6,1:3,i)*Rest'*l2d(:, i);                

    if (size(Rest, 1) > 0)
        Sigma2(1,1) = Sigma2(1,1) + Xsc(:, i)'*SigmasLines2D(:,:,i)*Xsc(:, i);
        Sigma2(2,2) = Sigma2(2,2) + Xec(:, i)'*SigmasLines2D(:,:,i)*Xec(:, i);
        Sigma2(1,2) = Sigma2(1,2) + Xsc(:, i)'*SigmasLines2D(:,:,i)*Xec(:, i);
        Sigma2(2,1) = Sigma2(2,1) + Xec(:, i)'*SigmasLines2D(:,:,i)*Xsc(:, i);
    end
    SigmaInv = inv(Sigma2);%Sigma2D_3D + 
    
    [Us,Ss,Vs] = svd(SigmaInv);
    SqrtSigmaInv = sqrt(Ss)*Vs';
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;    
    Mlines(row_ini:row_end,:) = SqrtSigmaInv*Mlines(row_ini:row_end,:);
end

%MtMe = Me'*Me;
MtMe = Mlines'*Mlines;
% Mpts = M'*M;

inds = [1 5 9 2 6 10 3 7 11 4 8 12];
% Me2 = Me1(inds, inds);
% indsRev = [1 4 7 10 2 5 8 11 3 6 9 12];
% Mpts = Mpts(indsRev, indsRev);
MtMe = MtMe + Me1;
% M = M'*M;
%     
[V1,S1]=eig(MtMe);

Km = V1(inds,4:-1:1);

% [V,S]=eig(M);
% 
% Km = V(:,4:-1:1);
end

