%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       In Submission to CVPR 2021
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

function [R, t] = optimize_noise_full(R0, t0, XX, xx, sigmas_2d, sigmas_3d, mode)
    if nargin < 7
        mode = 1;
    end    
    [rvec, jac] = rodrigues(R0);
    x = [rvec; t0];
%     fun = @(x) residual_full(x, XX, xx, sigmas_2d, sigmas_3d, R0, t0);
    n = size(sigmas_2d, 3);
    sigmas_2d_tr = zeros(n, 1);
    sigmas_3d_tr = zeros(n, 1);
    fun = @(x) residual_full(x, XX, xx, sigmas_2d, sigmas_3d);        
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'display', 'off');    
    x_fin = lsqnonlin(fun, x, zeros(0,0), zeros(0,0), OPTIONS);
    R = rodrigues(x_fin(1:3));
    t = x_fin(4:6);
end

function r = residual_full(x, XX, xx, sigmas_2d, sigmas_3d)
    R = rodrigues(x(1:3));
    t = x(4:6);
    n = size(XX, 2);
    XXc = R * XX + repmat(t, 1, n);
    xx_p = XXc(1:2, :) ./ repmat(XXc(3, :), 2, 1);
    dx_w = zeros(size(xx));
    dx = xx_p - xx;
    sigmas = zeros(n, 1);
    for i = 1:n
        jac = zeros(2, 3);
        jac(1,1) = 1.0 / XXc(3, i);
        jac(2, 2) = 1.0 / XXc(3, i);
        jac(1, 3) = -XXc(1, i) / XXc(3, i)^2;
        jac(2, 3) = -XXc(2, i) / XXc(3, i)^2;
        sigma = sigmas_2d(:, :, i) + jac * R * sigmas_3d(:, :, i) * R' * jac';                
        [U, S, V] = svd(sigma);        
        sigma_inv_sqrt = U * sqrt(inv(S)) * U';
        dx_w(:, i) = sigma_inv_sqrt * dx(:, i);
    end
    r = dx_w(:);        
end


function see = estimate_sigma_inv_sqrt(Xc, a, c)
    
    b = a / (1 + a * Xc(1:2)' * Xc(1:2));
    % b = 1 / (1/a + Xc(1:2)' * Xc(1:2))
    %1/a = c * Xc(3)^4 / sigma_p_2
    % b = 1 / (c * Xc(3)^4 / sigma_p_2 + Xc(1:2)' * Xc(1:2)) = 
    % = 1 / (sigma_d_2 / sigma_p_2 * Xc(3)^4 + \|Xc\|^2)
    alpha = 1.0 / norm(Xc(1:2)) * (1 - sqrt(1 - b));    
    see = 1.0/sqrt(c) * (eye(2) - alpha * Xc(1:2) * Xc(1:2)');
end