%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the CEPPnP 
% algorithms described in:
%
%       Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
%       Leveraging Feature Uncertainty in the PnP Problem. 
%       In Proceedings of BMVC, 2014. 
%
% Copyright (C) <2014>  <Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer>
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
% Luis Ferraz, CMTech-UPF, September 2014.
% luisferrazc@gmail.com,http://cmtech.upf.edu/user/62
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M, dM, Cw, Alph] = PrepareData(Pts,impts,Cw)

    if ~exist('Cw','var')
      Cw=define_control_points()';  
    end

    Xw=Pts';
    U=impts;

    %compute alphas (linear combination of the control points to represent the 3d points)
    Alph=compute_alphas(Xw,Cw');
    
    %Compute M
    M  =ComputeM(U(:),Alph);
    %Compute gradient of M
    dM =ComputedM(Alph);
end

function M = ComputeM(U,Alph)
    %ATTENTION U must be multiplied by K previously
    M = kron(Alph,[1 0 -1; 0 1 -1]);
    M(:,[[3,6,9,12]]) =  M(:,[3,6,9,12]) .* (U * ones(1,4));
    
end


function dM = ComputedM(Alph)
    dM = kron(Alph,[0 0 -1; 0 0 -1]);    
end