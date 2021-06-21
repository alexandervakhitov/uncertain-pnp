%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the REPPnP and EPPnP 
% algorithms described in:
%
%       Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
%       Very Fast Solution to the PnP Problem with Algebraic Outlier Rejection. 
%       In Proceedings of CVPR, 2014. 
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
% Luis Ferraz, CMTech-UPF, June 2014.
% luisferrazc@gmail.com,http://cmtech.upf.edu/user/62
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T,err] = EPPnP(Pts,impts)

    sol_iter = 1; %indicates if the initial solution must be optimized
    dims = 3;     %kernel dimensions
    
%    mPts = mean(Pts,2);
%    Pts  = Pts - (mPts * ones(1,size(Pts,2)));

    [M, Cw, Alph] = PrepareData(Pts,impts);
    
    t = sum((reshape(sum(M.^2),[3,4])));
    idx = find(t == min(t));
    %idx = find(sum((reshape(sum(M.^2),[3,4])))==0);
    if(~isempty(idx))
        M    = M(:,setdiff([1:12],[3*(idx-1)+1 3*(idx-1)+2 3*(idx-1)+3]));
        Cw   = Cw(:,setdiff([1:4],idx));
        Alph = Alph(:,setdiff([1:4],idx));
    end
    
    Km=kernel_noise(M,dims); %Compute kernel M
    
    [R, T, err] = KernelPnP(Cw, Km, dims, sol_iter);
%    T = T - R * mPts;

end

