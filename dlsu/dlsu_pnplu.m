function [R_est, t_est] = dlsu_pnplu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas = replace_with_trace(Sigmas3D);
    SigmasLines = zeros(6,6,0);% 
    l = zeros(3,0);
    Xs = zeros(3,0);
    Xe = zeros(3,0);
    xn = normalize_points(x, K);
    Sigmas2D = Sigmas2D / K(1,1)/K(1,1);    
    
    
    [R_est, t_est] = dlsu(X, xn, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, Sigmas2DLines, R_est_3, t_est_3);
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end
