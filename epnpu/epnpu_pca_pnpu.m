function [R_est, t_est] = epnpu_pca_pnpu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, mode)
    if nargin < 8
        mode = 2;
        m2 = 1;
    end
    [Sigmas3D, sigmas3d] = replace_with_trace(Sigmas3D);
    [Sigmas2D, sigmas2d] = replace_with_trace(Sigmas2D);
    n = size(X, 2);
    x3d_h = [X' ones(n,1)];
    x2d_h = [x' ones(n, 1)];    
%function [R,T,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,A,Sigmas3D,Sigmas2D,Rest,test,is_pca)    
    if mode == 2
        [R_est,t_est,Xc,best_solution,opt]=efficient_pnpu_gauss_fast(x3d_h,x2d_h,K,sigmas3d, sigmas2d, R_est_3, t_est_3, 0);
    else
        [R_est,t_est,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,K,Sigmas3D, Sigmas2D, R_est_3, t_est_3, 0);
    end
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
end