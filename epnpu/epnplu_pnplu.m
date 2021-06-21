function [R_est, t_est] = epnplu_pnplu(X, x, Sigmas3D, Sigmas2D, ...
    R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)
    sigmas2d = compute_traces(Sigmas2D);
    sigmas3d = compute_traces(Sigmas3D);
    sigmas3dlines = compute_traces(Sigmas3DLines);
%     Sigmas3D = replace_with_trace(Sigmas3D);
    n = size(X, 2);
%     x3d_h = [X' ones(n,1)];
%     x2d_h = [x' ones(n, 1)];    
%     xn = normalize_points(x, K);
    xn = x;
%EPnPLU(Xw,xx,xs, xe, Xs, Xe, Sigmas3D, Sigmas2D, l, Sigmas3DLines, Sigmas2DLines)
    [R_est, t_est, ff, s] = EPnPLU(X, xn, xs, xe, Xs, Xe, sigmas3d, ...
                        sigmas2d, l, sigmas3dlines, Sigmas2DLines, R_est_3, t_est_3);

    
%     [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end