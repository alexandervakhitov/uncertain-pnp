function [R_est, t_est] = run_p3p(XXw, xxn, f, R, t)
    p3p_inds = randsample(size(XXw,2),3);
    X_3 = XXw(:, p3p_inds);
    xu_3 = xxn(:, p3p_inds)/f;
    xuh_3 = [xu_3; ones(1,3)];        
    [R_est, t_est] = p3p_interface(X_3, xuh_3);   
    n_sol = size(t_est, 2);
    if n_sol == 0
        R_est = zeros(0, 0);
        t_est = zeros(0, 0);
        return;
    end
    pose_errors = zeros(n_sol,1);
    for sol_ind = 1:n_sol
        Rc = R_est(:,:,sol_ind);
        tc = t_est(:, sol_ind);
        pose_errors(sol_ind) = norm([Rc',  -Rc'*tc; 0 0 0 1]* [R t; 0 0 0 1] - eye(4));
    end        
    [best_err, best_sol] = min(pose_errors);
    R_est = R_est(:, :, best_sol);
    t_est = t_est(:, best_sol);
end