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

function [K, idinliers, i]=my_robust_kernel_noise(M,dimker, minerror)
 
    m   = size(M,1);
    id  =round(m/8);
    idx = 1:m;

    prev_sv = Inf;
    
    pairs = 1; %each correspondence is a couple of equations

    for i=1:10
    
        N = M(idx,:);
        [~,~,v] = svd(N'*N);
       
        if (pairs)
            error21    = M(1:2:end,:) * v(:,end);
            error22    = M(2:2:end,:) * v(:,end);
            error2     = sqrt(error21.^2 + error22.^2);
            
            [sv, tidx] = sort(error2);        

            med = sv(floor(m/8)); 

        else
            error2    = M * v(:,end);
            [sv, tidx] = sort(error2.^2);
            med = sv(floor(m/2)); 
        end
     
        ninliers = sum(sv<max(med,minerror));

        if (med >= prev_sv)
            break;
        else
            prev_sv = med;
            resv    = v;
            if(pairs)
                residx  = tidx(1:ninliers);
            else
                %always pairs = 1!! :P
               
            end
        end
        
        if(pairs)
            tidx2     = tidx'*2;
            ttidx     = [tidx2-1; tidx2];
            tidx2     = ttidx(:);
            idx       = tidx2(1:2*ninliers);
        else
            idx       = tidx(1:ninliers);
        end
    end
    
    K = resv(:,end-dimker+1:end);   
    idinliers = residx;
end