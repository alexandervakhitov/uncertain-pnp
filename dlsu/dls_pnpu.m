function [C_est_fin, t_est_fin] = dls_pnpu(X, x, K, Sigmas, Sigmas2D, R_est_3, t_est_3)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)    
    [C_est, t_est] = robust_dls_pnp(X, x);
    [C_est_fin, t_est_fin] = add_solution(C_est, t_est, R_est_3, t_est_3);
end