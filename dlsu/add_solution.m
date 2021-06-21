function [C_est_fin, t_est_fin] = add_solution(C_est, t_est, R_est_3, t_est_3)
    n_sol = size(t_est, 2);
    C_est_fin = zeros(3, 3, n_sol+1); 
    t_est_fin = zeros(3, n_sol+1);
    C_est_fin(:,:,1:n_sol) = C_est;
    C_est_fin(:,:,n_sol+1) = R_est_3;
    t_est_fin(:, 1:n_sol) = t_est;
    t_est_fin(:, n_sol+1) = t_est_3;
end