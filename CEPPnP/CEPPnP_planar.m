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

function [R,T, err] = CEPPnP(Pts,impts,Cu)
    
    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 3;     %kernel dimensions
   
    meanPts = mean(Pts,2);
    mPts = Pts-repmat(meanPts,1,size(Pts,2));
 
    [tu,~,~] = svd(mPts * mPts','econ');
    Cw = [tu, [0;0;0]];
    
    %Compute M
    [M, dM, Cw, Alph] = PrepareData2(Pts,impts, Cw);
    
    t = sum((reshape(sum(M.^2),[3,4])));
    idx = find(t == min(t));
    %idx = find(sum((reshape(sum(M.^2),[3,4])))<0.001);
    if(~isempty(idx))
        if (idx == 1)
            M    = M(:,[4:12]);
            dM   = dM(:,[4:12]);
            Cw   = Cw(:,[2,3,4]);
            Alph = Alph(:,[2,3,4]);
        elseif (idx ==2)
            M    = M(:,[1:3,7:12]);
            dM   = dM(:,[1:3,7:12]);
            Cw   = Cw(:,[1,3,4]);
            Alph = Alph(:,[1,3,4]);
        else
            M    = M(:,[1:6,10:12]);
            dM   = dM(:,[1:6,10:12]);
            Cw   = Cw(:,[1,2,4]);
            Alph = Alph(:,[1,2,4]);
        end
    end

     [~,~,Km] = svd(M'*M);
     [~,tKm,~]   = FNSani(M,dM,Cu,Km(:,end)); 
     Km = tKm(:,[end-dims+1:end]);
     [R, T, err] = KernelPnP(Cw, Km, dims, sol_iter);
     
     %to use the epnp method we need to project the points to the plane
     
  %   Pts2 = zeros(size(Pts));
  %   for i=1:size(Alph,1)
  %      Pts2(:,i) = sum(repmat(Alph(i,:),3,1) .* Cw,2);
  %   end
     
  %   err = 0;
  %   [R, T] = epnp_planar_solver(Pts2',impts',Alph, M, Km, Cw');
     
   % T = T - R * meanPts;
end

