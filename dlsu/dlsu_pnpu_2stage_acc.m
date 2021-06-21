function [R_est, t_est] = dlsu_pnpu_2stage_acc(X, x, K, Sigmas, Sigmas2D, R_est_0, t_est_0)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas_n = replace_with_trace(Sigmas);
    xn = normalize_points(x, K);
    Sigmas2D = Sigmas2D / K(1,1)/K(1,1);
    [R_est, t_est] = dlsu(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas_n, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_0, t_est_0);
    if (size(R_est, 3) == 1)
        [R_est_2, t_est_2] = dlsu(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est, t_est);
        if (size(R_est_2, 3) == 1)
            R_est = R_est_2;
            t_est = t_est_2;
        end
    end
    [R_est, t_est] = add_solution(R_est, t_est, R_est_0, t_est_0);
    
end