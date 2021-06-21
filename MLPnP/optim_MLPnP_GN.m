%     Steffen Urban email: urbste@googlemail.com
%     Copyright (C) 2016  Steffen Urban
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% 28.06.2016 Steffen Urban

function [Tout, statistic] = optim_MLPnP_GN(Tinit, points3D, ...
    rnull, snull, P, optimFlags)

% homogeneous to minimal
x = [Rodrigues2(Tinit(1:3,1:3))', Tinit(1:3,4)']';      

nrL = size(rnull,2);

% redundancy
redundanz = 2*nrL - length(x);    
% optim params
epsParam    = optimFlags.epsP;
epsFunc     = optimFlags.epsF;

% iteration params
cnt = 0;
stop = false;
invKll = P;
while cnt < optimFlags.maxit && stop == 0
    [r, J] = residualsAndJacobian(x, rnull, snull, points3D);                            
    % design matrix
    N = J.'*invKll*J;
    % System matrix
    g = J.'*invKll*r;

    dx = pinv(N)*g;
    if (max(abs(dx)) > 20 || min(abs(dx)) > 1)
        break;
    end
    dl = J*dx(1:end);

    if max(abs(dl)) < epsFunc || max(abs(dx(1:end))) < epsParam  
        x = x-dx;
        break;
    else
        % update parameter vector
        x = x-dx;              
    end  
    cnt = cnt+1;
end % while loop


% minimal to homogeneous
Tout = [Rodrigues2(x(1:3)) x(4:6)];

% empirical variance factor
resV = r.'*invKll*r;
if redundanz > 0
    if redundanz < nrL
        s0 = 1;  
    else
        s0 = resV / redundanz;
    end
else
    s0 = NaN;
end
% variance-covariance matrix
Qxx = pinv(N);       
% cofactor matrix of "adjusted observations"
Qldld = J*Qxx*J';

statistic = {resV, r, J, Qxx, s0, Qldld, sqrt(s0.*diag(Qxx))};
end