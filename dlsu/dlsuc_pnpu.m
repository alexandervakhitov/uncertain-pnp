function [R_est, t_est] = dlsuc_pnpu(X, x, K, Sigmas, Sigmas2D, R_est_3, t_est_3)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    xn = normalize_points(x, K);
    Sigmas2D = Sigmas2D / K(1,1)/K(1,1);
    [R_est, t_est] = dlsuc(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3, t_est_3);
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end