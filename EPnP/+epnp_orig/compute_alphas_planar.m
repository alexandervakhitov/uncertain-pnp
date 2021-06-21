function [Alph Alph_c ma] = compute_alphas_planar(Xw,Cw)

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


n = size(Xw, 2); %number of 3d points

%generate auxiliar matrix to compute alphas
C = [Cw(:, 1:2)'; ones(1,3)];
X=[Xw(1:2, :); ones(1,n)];
Alph = inv(C)*X;
ma = mean(Alph, 2);
Alph_c = Alph - repmat(ma, 1, n);


