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
% This script contains code for the DLS(L)U method
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


function [R_est, t_est] = dlsu_full(X, x, Sigmas, Sigmas2D, R_est_3, t_est_3, ...
                                    Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, ...
                                    full_cov_mode)

    if nargin < 12
        full_cov_mode = 0;
    end
    
    if length(size(Sigmas)) < 3
        % isotropic covariances of type s^2 I, only s^2 are provided
        npt = size(X, 2);
        Sigmas = repmat(reshape(eye(3), 3, 3, 1), 1, 1, npt);
        Sigmas2D = repmat(reshape(eye(2), 2, 2, 1), 1, 1, npt);
        t_est_3 = 1;
    end
    
    %centralize the 3D points
    meanX = mean(X,2);
    X = X - repmat(meanX, 1, size(X,2));    
    if length(t_est_3) > 1
        t_est_3 = t_est_3 - R_est_3 * meanX;
    end
    
    if full_cov_mode == 0
        %replace covariance with isotropic aproximations
        sigmas2 = compute_traces(Sigmas);
        sigmas2d_2 = compute_traces(Sigmas2D);        
        sigmas_lines_2 = compute_traces(Sigmas3DLines);            
    end
    if (size(Xs, 1) == 0)
        Xs = zeros(3, 0);
        Xe = zeros(3, 0);
        Sigmas3DLines = zeros(6, 6, 0);
    else
        %centralize 3D lines
        Xs = Xs - repmat(meanX, 1, size(Xs, 2));
        Xe = Xe - repmat(meanX, 1, size(Xe, 2));
    end
    tic;
    xn = x;       
    if full_cov_mode == 0
        [R_est, t_est, costs_est] = dlsu_fast(X, xn, Xs, Xe, l, sigmas2, ...
                                              sigmas_lines_2, sigmas2d_2, ...
                                              Sigmas2DLines, R_est_3, t_est_3);
    else
        [R_est, t_est, costs_est] = dlsu(X, xn, Xs, Xe, l, Sigmas, ...
                                         Sigmas3DLines, Sigmas2D, ...
                                         Sigmas2DLines, R_est_3, t_est_3);
    end
    %choose solutions with low cost value
    m = costs_est < 0.1;
    R_est = R_est(:, :, m);
    t_est = t_est(:, m);

    if (size(R_est, 3) == 0)
        % try random rotations of 3D features
        cnt = 0;
        R_est_r = zeros(3, 3, 100);
        t_est_r = zeros(3, 100);
        costs_r = zeros(100, 1);
        ind = 1;
        while cnt < 3
            randax = 3*randn(3,1);
            randax(randi(3)) = 0;
            Rr = rodrigues(randax);
            if size(R_est_3, 1) > 0
                R_est_new = R_est_3*Rr';
            else
                R_est_new = R_est_3;            
            end
            if size(Xs, 1) == 0
                Xs = zeros(3, 0);
                Xe = zeros(3, 0);
            end
            t_est_new = t_est_3;
            if full_cov_mode == 0
                [R_est, t_est, costs_est] = dlsu_fast(Rr*X, xn, Rr*Xs, ...
                                                      Rr*Xe, l, sigmas2, ...
                                                      sigmas_lines_2, ...
                                                      sigmas2d_2, ...
                                                      Sigmas2DLines, ...
                                                      R_est_new, t_est_new);
            else
                n = size(Sigmas, 3);
                SigmasR = zeros(3, 3, n);
                Sigmas3DLinesR = zeros(6, 6, n);
                if full_cov_mode == 1
                    for j = 1:n
                        SigmasR(:, :, j) = Rr * Sigmas(:, :, j) * Rr';
                    end
                    n_ln = size(Sigmas3DLines, 3);
                    for j = 1:n_ln
                        Sigmas3DLinesR(1:3, 1:3, j) = Rr * Sigmas3DLines(1:3, 1:3, j) * Rr';
                        Sigmas3DLinesR(4:6, 4:6, j) = Rr * Sigmas3DLines(4:6, 4:6, j) * Rr';
                    end                
                else
                    SigmasR = Sigmas;
                    Sigmas3DLinesR = Sigmas3DLines;
                end
                [R_est, t_est, costs_est] = dlsu(Rr*X, xn, Rr*Xs, Rr*Xe, l, ...
                                                 SigmasR, Sigmas3DLinesR, ...
                                                 Sigmas2D, Sigmas2DLines, ...
                                                 R_est_new, t_est_new);
            end
            for k = 1:size(R_est,3)
                R_est_r(:, :, ind) = R_est(:,:,k)*Rr;
                t_est_r(:, ind) = t_est(:,k);
                costs_r(ind) = costs_est(k);
                ind = ind+1;
            end 
            cnt = cnt+1;
        end
        [~, min_ind] = min(costs_r(1:ind-1));
        R_est = R_est_r(:, :, 1:ind-1);
        t_est = t_est_r(:, 1:ind-1);
    end   
        
    for si = 1:size(R_est, 3)
        t_est(:, si) = t_est(:, si) - R_est(:, :, si) * meanX;                    
    end
            
end