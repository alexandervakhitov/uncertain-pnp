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

function [R, t] = optimize_noise_2d(R0, t0, XX, xx, sigmas_2d, Xs, Xe, ll, sigmas_lines_2d)
    if (nargin < 6)
        Xs = [];
        Xe = [];
        ll = [];
        sigmas_lines_2d = [];
    end
    if size(Xs, 1) == 0
        Xs = zeros(3, 0);
        Xe = zeros(3, 0);
        ll = zeros(3, 0);
    end
    [rvec, jac] = rodrigues(R0);
    x = [rvec; t0];
    sigmas_2d_inv_sqrt = zeros(size(sigmas_2d));
    n = size(XX, 2);
    for i = 1:n
%         [U, S, V] = svd(sigmas_2d(:, :, i));                
%         sigmas_2d_inv_sqrt(:, :, i) = U * sqrt(inv(S)) * U';
        sigma_inv = inv(sigmas_2d(:, :, i));
        sigmas_2d_inv_sqrt(:, :, i) = chol(sigma_inv);
    end    
    
    sigmas_2d_lines_inv_sqrt = sqrt(1./sigmas_lines_2d);
    
    fun = @(x) residual_2d(x, XX, xx, Xs, Xe, ll, sigmas_2d_inv_sqrt, sigmas_2d_lines_inv_sqrt);
    OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'display', 'off', 'SpecifyObjectiveGradient', true);        %'SpecifyObjectiveGradient',true 'StepTolerance', 1e-9, 'FunctionTolerance', 1e-9, 
    x_fin = lsqnonlin(fun, x, zeros(0,0), zeros(0,0), OPTIONS);    
    R = rodrigues(x_fin(1:3));
    t = x_fin(4:6);
end

