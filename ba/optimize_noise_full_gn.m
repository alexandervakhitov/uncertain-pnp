%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       A. Vakhitov, L. Ferraz, A. Agudo, F. Moreno-Noguer
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       CVPR 2021
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, t, statistic] = optimize_noise_full_gn(R0, t0, XX, xx, sigmas_2d, sigmas_3d, ...
    Xs, Xe, ll, sigmas_lines_2d, sigmas_lines_3d, optimFlags)

    if size(Xs, 1) == 0
       Xs = zeros(3, 0);
       Xe = zeros(3, 0);
       ll = zeros(3, 0);
    end
    [rvec, jac] = rodrigues(R0);
    x = [rvec; t0];    
    n = size(XX, 2);
    nl = size(Xs, 2);
    nf = n + nl;
    sigmas_2d_inv = zeros(2*nf, 2*nf);
    
    % redundancy
    redundanz = 2*n - length(x);    
    % optim params
    epsParam    = optimFlags.epsP;
    epsFunc     = optimFlags.epsF;

    % iteration params
    cnt = 0;
    stop = false;
    invKll = sigmas_2d_inv;
    sigma_proj_3d = zeros(2, 2, nl);
    while cnt < optimFlags.maxit && stop == 0
        [r, J] = residual_2d(x, XX, xx, Xs, Xe, ll);                          
        
        R = rodrigues(x(1:3));
        for i = 1:n        
            invKll(2*i-1:2*i, 2*i-1:2*i) = inv(sigmas_2d(:, :, i) + J(2*i-1:2*i, 4:6) * R * sigmas_3d(:, :, i) * R' * J(2*i-1:2*i, 4:6)');
        end    
        for i = n+1:n+nl
            li = i - n;                        
            sigma_proj_3d(1, 1, li) = J(2*i-1, 4:6) * R * sigmas_lines_3d(1:3, 1:3, li) * R' * J(2*i-1, 4:6)';
            sigma_proj_3d(1, 2, li) = J(2*i-1, 4:6) * R * sigmas_lines_3d(1:3, 4:6, li) * R' * J(2*i, 4:6)';
            sigma_proj_3d(2, 1, li) = sigma_proj_3d(1, 2);
            sigma_proj_3d(2, 2, li) = J(2*i, 4:6) * R * sigmas_lines_3d(4:6, 4:6, li) * R' * J(2*i, 4:6)';            
            invKll(2*i-1:2*i, 2*i-1:2*i) = inv(sigmas_lines_2d(li) * eye(2) + sigma_proj_3d(:, :, li));
        end
        
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
    R = rodrigues(x(1:3));
    t = x(4:6);

    % empirical variance factor
    resV = r.'*invKll*r;
    if redundanz > 0
        if redundanz < n
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

