function [R_est, t_est] = epnpl_pnplu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)    
    n = size(X, 2);
    xn = normalize_points(x, K);
%function [R,t,ff, s] = EPnPLS_GN(Xw,xx,xs, xe, Xs, Xe, l)
%     l = zeros(3,0);
%     xs = zeros(2,0);
%     xe = zeros(2,0);
%     Xs = zeros(3,0);
%     Xe = zeros(3,0);
    
    [R_est, t_est, ff, s] = EPnPLS_GN(X, xn, xs, xe, Xs, Xe, l);

    
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end