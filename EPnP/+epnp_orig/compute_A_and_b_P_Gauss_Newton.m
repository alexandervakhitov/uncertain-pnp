function [A,b]=compute_A_and_b_P_Gauss_Newton(current_betas,rho,L)

% COMPUTE_A_AND_B_GAUSS_NEWTON  
%
%       A: Jacobian
%       b: error vector
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



A=zeros(3,3);
b=zeros(3,1);

cb=current_betas;

B=[cb(1)*cb(1),
   cb(1)*cb(2),
   cb(2)*cb(2),
   cb(1)*cb(3),
   cb(2)*cb(3),
   cb(3)*cb(3)];

for i=1:3
    A(i,1)=2*cb(1)*L(i,1)+cb(2)*L(i,2)+cb(3)*L(i,4);
    A(i,2)=cb(1)*L(i,2)+2*cb(2)*L(i,3)+cb(3)*L(i,5);
    A(i,3)=cb(1)*L(i,4)+cb(2)*L(i,5)+2*cb(3)*L(i,6);            
    b(i)=(rho(i)-L(i,:)*B);    
end
    



    

