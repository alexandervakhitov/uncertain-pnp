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
% AV added pre/post processing with centering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T, err] = CEPPnP(Pts_init,impts,Cu)

    Pts = Pts_init - repmat(mean(Pts_init,2),1,size(Pts_init,2));    
  
    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 4;     %kernel dimensions
   
    %Compute M
    [M, dM, Cw, Alph] = PrepareData2(Pts,impts);
  
     [~,~,Km] = svd(M'*M);
     [~,tKm,~]   = FNSani(M,dM,Cu,Km(:,end)); 
     Km = tKm(:,[end-dims+1:end]);
     [R, T, err] = KernelPnP(Cw, Km, dims, sol_iter); 
     
     T = T - R * mean(Pts_init,2);
end

