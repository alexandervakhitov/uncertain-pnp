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
% This script contains the 2D noise experiments 
% from the main paper (Figure 3, top row).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
IniToolbox;

% camera parameters
width = 640;
height = 480;
f= 800;
K = [f 0 0
     0 f 0
     0 0 1];

% experimental parameters
nl = [1,2,3,4,5,6,7,8,9,10];

%percentage of samples for each sigma
nlsamples = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; 
npts = [10:10:110];
num = 40;

% compared methods
A = zeros(size(npts));
B = zeros(num,1);

name   = {'EPnP+GN', 'DLS',           'OPnP', 'CEPPnP', 'MLPnP'  ,       'DLSU',       'DLSU*',   'EPnPU',     'EPnPU*'};
fun    = {@EPnP_GN,  @robust_dls_pnp, @OPnP,  @CEPPnP,  @MLPNP_with_COV, @dlsu_full, @dlsu_full,  @epnpu_full, @epnpu_full};
marker = {'-.',      '-.',            '-.',   '--',     '-',             '-',          '-.',      '-',        '-.' };
color  = {'c',       'm',             'b'     'k',      'k',             'g',          'g',       'r',         'r' };
 
method_list= struct('name', name, 'f', fun, 'mean_r', A, 'mean_t', A, ...
                    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, ...
                    'r', B, 't', B, 'marker', marker, 'color', color, ...
                    'markerfacecolor', color);

% experiments
for i=1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d (num sg = %d ): ', npt, length(nl));
    
     for k= 1:length(method_list)
        method_list(k).c = zeros(1,num);
        method_list(k).e = zeros(1,num);
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    index_fail = cell(1,length(name));    
    for j= 1:num
        
        % generate 3d coordinates in camera space
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
        avg_depth = 6;
%         Xc= [xrand(1,npt,[-8 8]); xrand(1,npt,[-8 8]); xrand(1,npt,[16 32])];
%         avg_depth = 24;

        % generate camera
        t = mean(Xc,2);
        R = rodrigues(randn(3,1));
        XXw = inv(R)*(Xc-repmat(t,1,npt));        
        
        % project points
        [xxn, v, Cu, cov] = project_points_compute_covs(Xc, f, nl, nlsamples);
        
        % 3D covariances are zero
        Sigmas3D = zeros(3,3,npt);
                     
        %call p3p to get R_est, t_est
        [R_est, t_est] = run_p3p(XXw, xxn, f, R, t);
        if size(R_est, 1) < 3
            continue;
        end
        
        % pose estimation
        for k= 1:length(method_list)            
            [is_fail, R1, t1, tcost] = run_pnpl_method(method_list(k), XXw, xxn, v, ...
                                                       f, Cu, cov, Sigmas3D, ...
                                                       avg_depth, R_est, t_est);
            if (is_fail)
                index_fail{k} = [index_fail{k}, j];
                continue;
            end
            %choose the solution with the smallest error 
            [error, y, best_id, ercorr] = choose_solution(R1, t1, R, t, XXw, Xc);
            
            method_list(k).c(j)= tcost * 1000;
            method_list(k).e(j)= ercorr;
            method_list(k).r(j)= y(1);
            method_list(k).t(j)= y(2);
                        
        end

        showpercent(j,num);
    end
    fprintf('\n');
    
    method_list = save_experiment(method_list, i, index_fail, num);
end

mkdir('results/');
exp_lbl = 'results/ordinary2Dresults';
save(exp_lbl, 'method_list', 'npts');
save_pref = 'ord2d';
plotOrdinary3D_2Dnoise(exp_lbl, save_pref);