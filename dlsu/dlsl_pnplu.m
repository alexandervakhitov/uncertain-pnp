function [R_est, t_est] = dlsl_pnplu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)
    xn = normalize_points(x, K);    
%     [Rs, ts] = dlsl(X, x, Xs, Xe, l)
    X = zeros(3,0);
    xn = zeros(2,0);
    [R_est, t_est] = dlsl(X, xn, Xs, Xe, l);
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end