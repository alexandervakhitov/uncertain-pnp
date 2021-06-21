function Me1 = compute_M1_planar(U,Alph,A)

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

At = Alph;
Ut = U;
M1 = [fu*ones(n, 3) zeros(n, 3) repmat(u0-Ut(:, 1), 1, 3)].*repmat(At, 1, 3);
M2 = [zeros(n, 3) fv*ones(n, 3) repmat(v0-Ut(:, 2), 1, 3)].*repmat(At, 1, 3);

Me1 = M1'*M1 + M2'*M2;


