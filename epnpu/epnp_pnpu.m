function [R_est, t_est] = epnp_pnpu(X, x, K, Sigmas, Sigmas2D, R_est_3, t_est_3)
    n = size(X, 2);
    x3d_h = [X' ones(n,1)];
    x2d_h = [x' ones(n, 1)];    
    [R_est,t_est,Xc,best_solution,opt]=efficient_pnp_gauss(x3d_h,x2d_h,K,false);
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
end