%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the CEPPnP 
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

clear; clc;
IniToolbox;

% camera's parameters
width= 640;
height= 480;
f= 800;
K = [f 0 0
     0 f 0
     0 0 1];


% experimental parameters
nl_2d = [1,2,3,4,5,6,7,8,9,10];
nl = [1,2,3,4,5,6,7,8,9,10]/40;

%nl= 5*ones(1,10);
nlsamples = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %percentatge of samples for each sigma
npts = [30:10:110];
isotropic = false;
num = 400;

% compared methods
A= zeros(size(npts));
B= zeros(num,1);

funs = {[],    @EPnP_GN,  @MLPNP_with_COV, @epnpu_full};
lbls = {'P3P', 'EPnP', 'MLPnP', 'EPnPU'};
colors = {'m', 'c',    'k',   'r'};

markers = {'-.', ':'};
mod_lbl = {'URef', 'SRef'};
rtypes = {1, 1};
ref_types = {1, 3};

ind = 1;
for meti = 1:length(funs)
    for modi = 1:length(mod_lbl)
        name{ind} = lbls{meti};        
        name{ind} = [name{ind} '+' mod_lbl{modi}];        
        fun{ind} = funs{meti};
        marker{ind} = markers{modi};
        color{ind} = colors{meti};
        markerfacecolor{ind} = colors{meti};
        ransac_type{ind} = rtypes{modi};
        ref_type{ind} = ref_types{modi};
        ind = ind + 1;
    end
end

method_list= struct('name', name, 'f', fun, 'mean_r', A, 'mean_t', A,...
    'med_r', A, 'med_t', A, 'std_r', A, 'std_t', A, 'r', B, 't', B,...
    'marker', marker, 'color', color, 'markerfacecolor', markerfacecolor, ...
    'ransac_type', ransac_type, 'refinement_type', ref_type);

% experiments
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d (num sg = %d ): ', npt, length(nl));
   
    
     for k= 1:length(method_list)
        method_list(k).c = zeros(1,num);
        method_list(k).e = zeros(1,num);
        method_list(k).r = zeros(1,num);
        method_list(k).t = zeros(1,num);
    end
    
    %index_fail = [];
    index_fail = cell(1,length(name));
    p3p_fail_cnt = 0;
    
    for j= 1:num
        
        % generate 3d coordinates in camera space
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
        avg_depth = mean(Xc(3, :));
            
        t= mean(Xc,2);
        R= rodrigues(randn(3,1));
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        XXwt = XXw;
        
        % generate 3D points
        [XXw, Sigmas3D] = generate_3d_points_with_covs(XXw, nl, nlsamples, isotropic);
        
        % project points
        [xxn, v, Cu, cov, nls2d] = project_points_compute_covs(Xc, f, nl_2d, nlsamples, true);         
               
        % RANSAC
        timer_ransac = tic;
        fittingfn = @p3p_interface;
        n_sample = 3;%sample size
        scales = nls2d;
        thr = 6;
        if size(xxn, 2) < n_sample                
            index_fail{k} = [index_fail{k}, j];
            continue;
        end
        
        [R_est, t_est, inliers, trialcount] = ransac(K, XXw, xxn, fittingfn, n_sample, thr, scales);        
        
        if method_list(k).ransac_type == 2
            thr_3d = 25;
            [inliers, R_est, t_est] = ransac_eval_3d(K, R_est, t_est, XXw, xxn, thr, ...
                    thr_3d, scales, Sigmas3D);                
        end
        
        cost_ransac = toc(timer_ransac);

        if size(R_est, 1) == 0 || length(inliers) < 6
            %solution not found; go for the next trial
                index_fail{k} = [index_fail{k}, j];
            continue;
        end

        XXwi = XXw(:, inliers);
        xxni = xxn(:, inliers);
        Cui = Cu(:, :, inliers);
        covi = cov(:, inliers);
        vi = v(:, inliers);
        Sigmas_3di = Sigmas3D(:, :, inliers);
        
        % pose estimation
        for k= 1:length(method_list)            

            [is_fail, R1, t1, tcost] = run_pnpl_method(method_list(k), XXwi, xxni, vi, ...
                                                       f, Cui, covi, Sigmas_3di, ...
                                                       avg_depth, R_est, t_est);
            if (is_fail)
                index_fail{k} = [index_fail{k}, j];
                continue;
            end
            
            [inliers_n, R0, t0] = ransac_eval(K, R1, t1, XXw, xxn, thr, scales);
            
            if method_list(k).ransac_type == 2
                thr_3d = 25;
                [inliers_n, R0, t0] = ransac_eval_3d(K, R1, t1, XXw, xxn, thr, ...
                        thr_3d, scales, Sigmas3D);                
            end
            if length(inliers_n) < 3
                index_fail{k} = [index_fail{k}, j];
                continue;
            end
            inliers = inliers_n;

            timer_ref = tic;
            if method_list(k).refinement_type == 1

                XXwf = XXw(:, inliers);
                xxnf = xxn(:, inliers) / f;
                Cunf = Cu(:, :, inliers);
                Sigmas3Df = Sigmas3D(:, :, inliers);
                optimFlags.epsP  = 1e-6;
                optimFlags.epsF  = 1e-6;
                optimFlags.maxit = 5;
                optimFlags.tau   = 1e-4;
                [Rf, tf, statistic] = optimize_noise_full_gn(R0, t0, XXwf, xxnf, ...
                    Cunf, Sigmas3Df, [], [], [], [], [], optimFlags);
            elseif method_list(k).refinement_type == 2
                XXwf = XXw(:, inliers);
                vf = v(:, inliers);
                covf = cov(:, inliers);
                [Rf, tf] = MLPNP_refine(R0, t0, XXwf, vf, covf);                                
            elseif method_list(k).refinement_type == 4
                XXwf = XXw(:, inliers);
                vf = v(:, inliers);
                covf = cov(:, inliers);
                Rf = R0;
                tf = t0;
            elseif method_list(k).refinement_type == 5
                XXwf = XXw(:, inliers);
                xxnnf = xxn(:, inliers) / f;
                Cuf = Cu(:, :, inliers);
                Cunf = Cu_new(:, :, inliers);   
                Sigmas3Df = Sigmas3D(:, :, inliers);
                [Rf, tf] = optimize_noise_full(R0, t0, XXwf, xxnnf, Cunf, Sigmas3Df);                
            else
                XXwf = XXw(:, inliers);
                xxnf = xxn(:, inliers)/f;
                Cunf = Cu(:, :, inliers);                                
                [Rf, tf] = optimize_noise_2d(R0, t0, XXwf, xxnf, Cunf);                
            end
            cost_nl = toc(timer_ref);
            R1 = zeros(3, 3, 1);
            R1(:,:,1) = Rf;
            t1 = zeros(3, 1);
            t1(:, 1) = tf;                    
            
            %choose the solution with smallest error 
            [error, y, best_id, ercorr] = choose_solution(R1, t1, R, t, XXw, Xc);
            cost  = tcost +  cost_nl;
            method_list(k).c(j)= cost * 1000;
            method_list(k).e(j)= ercorr;
            method_list(k).r(j)= y(1);
            method_list(k).t(j)= y(2);            
                                        

        end

        showpercent(j,num);
    end
    fprintf('\n');
    
    % save result
    method_list = save_experiment(method_list, i, index_fail, num);
end

save results/ordinary3DresultsRef method_list npts;

plotOrdinary3D('results/ordinary3DresultsRef', 'ord3d_ref_', true);