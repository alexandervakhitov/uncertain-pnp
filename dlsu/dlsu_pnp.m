function poses = dlsu_pnp(X, x, Sigmas, Sigmas2D, R_est, t_est)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    poses = dlsu(X, x, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est, t_est);
end