function [R_est, t_est] = dls_pnplu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas = replace_with_trace(Sigmas3D);
    SigmasLines = replace_with_trace(Sigmas3DLines);
    xn = normalize_points(x, K);
%     ln = normalize_lines(l, K);
    Sigmas2D = Sigmas2D / K(1,1)/K(1,1);        
    %dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Xs = zeros(3,0);
    Xe = zeros(3,0);
    n = size(X, 2);
    Sigmas = zeros(3,3,n);
    for i = 1:n
        Sigmas2D(:,:,i) = eye(2);
    end
    [R_est, t_est] = dlsu(X, xn, Xs, Xe, zeros(3,0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3, t_est_3);
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end