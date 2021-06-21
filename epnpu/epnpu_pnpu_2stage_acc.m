function [R_est, t_est] = epnpu_pnpu_2stage_acc(X, x, K, Sigmas3D, Sigmas2D, R_est_0, t_est_0)
    Sigmas3D_n = replace_with_trace(Sigmas3D);
    n = size(X, 2);
    x3d_h = [X' ones(n,1)];
    x2d_h = [x' ones(n, 1)];    
%function [R,T,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,A,Sigmas3D,Sigmas2D,Rest,test,is_pca)    
    [R_est,t_est,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,K,Sigmas3D_n, Sigmas2D, R_est_0, t_est_0, true);
    if (size(R_est, 3) == 1)
        [R_est_2,t_est_2,Xc,best_solution,opt]=efficient_pnpu_gauss(x3d_h,x2d_h,K,Sigmas3D, Sigmas2D, R_est, t_est, true);
        if (size(R_est_2,3) == 1)
            R_est = R_est_2;
            t_est = t_est_2;
        end
    end
    [R_est, t_est] = add_solution(R_est, t_est, R_est_0, t_est_0);
end