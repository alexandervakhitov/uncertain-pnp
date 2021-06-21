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
%
% This script contains code for the EPnP(L)U method
%
% X - 3 x n_pt array of 3D points
% x - 2 x n_pt array of 2D point projections
% Sigmas - 3 x 3 x n_pt array of 3D point covariances
% Sigmas2D - 2 x 2 x n_pt array of 2D point covariances
% R_est_3 - 3 x 3 rotation matrix of a pose hypothesis, or a 0x0 matrix 
%           if none is provided
% t_est_3 - 3 x 1 translation vector of a pose hypothesis, or a single
%           float value of the average scene depth 
% Xs, Xe - 3 x n_l arrays of 3D line endpoints
% l - 3 x n_l array of line equations in the image plane, 
%     l_i = [c_i, s_i, d], c_i^2 + s_i^2 = 1
% Sigmas3DLines - 6 x 6 x n_l array of 3D covariances of the line endpoints
% Sigmas2DLines - n_l array of 2D covariances of the line endpoints
% full_cov_mode - whether we use a method with a pose hypothesis (full_cov_mode = 1),
%                 or not (full_cov_mode = 0)
%

function [R_est, t_est] = epnpu_full(X, x, Sigmas3D, Sigmas2D, ...
                                     R_est_3, t_est_3, Xs, Xe, l, ...
                                     Sigmas3DLines, Sigmas2DLines, ...
                                     full_cov_mode)    
    
    K = eye(3);
    if nargin < 5
        R_est_3 = [];
        t_est_3 = [];
    end
    if nargin < 12
        full_cov_mode = 0;
    end

    if (size(Xs, 2) == 0)
        Xs = zeros(3, 0);
        Xe = zeros(3, 0);
        l = zeros(3, 0);
        Sigmas3DLines = zeros(6, 6, 0);
        Sigmas2DLines = zeros(0, 0);
    end
    tic;                
    
    n = size(X, 2);
    x3d_h = [X' ones(n,1)];
    x2d_h = [x' ones(n, 1)];    
    
    if full_cov_mode == 0        
        sigmas3d = compute_traces(Sigmas3D);        
        sigmas2d = compute_traces(Sigmas2D);        
        sigmas3dlines = compute_traces(Sigmas3DLines);
        if size(Xs, 1) == 0
            Xs = zeros(3, 0);
            Xe = zeros(3, 0);
        end
        [R_est,t_est,~,~,~] = efficient_pnpu_gauss_fast(x3d_h,x2d_h,K,...
                                                        sigmas3d, sigmas2d, ...
                                                        Xs, Xe, l, ...
                                                        sigmas3dlines, Sigmas2DLines, ...
                                                        R_est_3, t_est_3);
    else

        [R_est,t_est,~,~,~,~] = efficient_pnpu_gauss_2(x3d_h,x2d_h, K, ...
                                                       Sigmas3D, Sigmas2D, ...
                                                       Xs, Xe, l, Sigmas3DLines, ...
                                                       Sigmas2DLines, ...
                                                       R_est_3, t_est_3);

    end
end