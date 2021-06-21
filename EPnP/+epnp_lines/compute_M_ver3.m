function Km = compute_M_ver3(U, Alph,A, Alph_S, Alph_E, l2d, wl)

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

Me1 = M1'*M1 + M2'*M2;
%from here works only for u0 = v0 = 0, fu = fv

AS = [Alph_S Alph_E]';
l22 = repmat(l2d', 2, 1);
Me = [repmat(l22(:,1),1,4) repmat(l22(:,2),1,4) repmat(l22(:,3),1,4) ] .* repmat(AS, 1, 3);


MtMe = Me'*Me;
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

