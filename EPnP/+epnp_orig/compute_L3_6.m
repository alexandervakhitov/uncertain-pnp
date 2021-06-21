function L=compute_L3_6(K)

% COMPUTE_L3_6  
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

L=zeros(3,6);

ind = 1;
for i = 1:3
    for j = i+1:3
        P1 = K(3*i-2:3*i, :);
        P2 = K(3*j-2:3*j, :);
        D = P1-P2;
        DTD = D'*D;
        %b1^2 b1b2 b2^2 b1b3 b2b3 b3^2
        L(ind, 1) = DTD(1,1);
        L(ind, 2) = 2*DTD(1,2);
        L(ind, 3) = DTD(2,2);
        L(ind, 4) = 2*DTD(1,3);
        L(ind, 5) = 2*DTD(2, 3);
        L(ind, 6) = DTD(3,3);   
        ind = ind+1;
    end
end
