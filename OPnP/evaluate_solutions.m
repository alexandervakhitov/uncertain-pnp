function [ err err_list] = evaluate_solutions( sols,groundtruth )
% evalute solutions by calculating l2 norm between solutions and
% groundtruth
sols = [sols -sols];
sols_dim = size(sols,2);

err_list  = ( sqrt (abs ( sum( (real(sols) - repmat(groundtruth,1,sols_dim) ) .^2 ) ) ) ) /norm(groundtruth) ;  

err  = min(err_list);
