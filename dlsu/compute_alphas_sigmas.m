function [Alph, Sigmas3DAlph] =compute_alphas_sigmas(Xw,Sigmas3D,Cw)

% COMPUTE_ALPHAS Barycentric coordinates computation
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
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 


n=size(Xw,1); %number of 3d points

%generate auxiliar matrix to compute alphas
C=[Cw';ones(1,4)];
X=[Xw';ones(1,n)];
CI = inv(C);
Alph_=CI*X;
n = size(Sigmas3D,3);
CI = inv(C);
Sigmas3DAlpha = zeros(3, 3, n);
for i = 1:n
    Sigmas3DAlpha(:,:,i) = CI*Sigmas3D(:,:,i)*CI';
end

Alph=Alph_';
