function [R_est, t_est] = epnpu_pnpu(X, x, Sigmas3D, Sigmas2D, R_est_3, t_est_3, pca_mode)
    K = eye(3);
    mode = 2;
    if nargin < 5
        R_est_3 = [];
        t_est_3 = [];
    end
    if nargin < 7
        pca_mode = 1;
    end
    tic;                
    
    tc1 = toc;
    n = size(X, 2);
    x3d_h = [X' ones(n,1)];
    x2d_h = [x' ones(n, 1)];    
%function [R,T,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,A,Sigmas3D,Sigmas2D,Rest,test,is_pca)    
    tic;
    if mode == 2        
        sigmas3d = compute_traces(Sigmas3D);
        sigmas2d = compute_traces(Sigmas2D);
        [R_est,t_est,Xc,best_solution,opt]=efficient_pnpu_gauss_fast(x3d_h,x2d_h,K,sigmas3d, sigmas2d, R_est_3, t_est_3, pca_mode);
    else
        [Sigmas3D, sigmas3d] = replace_with_trace(Sigmas3D);
        [Sigmas2D, sigmas2d] = replace_with_trace(Sigmas2D);
        [R_est,t_est,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,K,Sigmas3D, Sigmas2D, R_est_3, t_est_3, true);
    end
%     [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    tc2 = toc;
%      fprintf('%f / %f\n', tc1, tc2);
end