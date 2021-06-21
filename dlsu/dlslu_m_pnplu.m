function [R_est, t_est] = dlslu_m_pnplu(X, x, Sigmas, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, mode)
% function [R_est, t_est] = dlsu_m_pnpu(X, x, K, Sigmas, Sigmas2D, R_est_3, t_est_3, mode)
    if nargin < 13
        mode = 1;
    end
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas0 = Sigmas;
    Sigmas2D_0 = Sigmas2D;
    [Sigmas, sigmas2] = replace_with_trace(Sigmas);
    [Sigmas2D, sigmas2d_2] = replace_with_trace(Sigmas2D);
    sigmas_lines_2 = compute_traces(Sigmas3DLines);        
    n = size(X, 2);    
    
    xn = x;
    [R_est, t_est] = dlsu_fast(X, xn, Xs, Xe, l, sigmas2, sigmas_lines_2 , sigmas2d_2, Sigmas2DLines, R_est_3, t_est_3);
    if (size(R_est, 3) == 0)
        cnt = 0;
        while cnt < 5 && size(R_est, 3) == 0
            Rr = rodrigues(3*randn(3,1));        
            
            %[R_est, t_est] = dlsu(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3*Rr', t_est_3);
            [R_est, t_est] = dlsu_fast(Rr*X, xn, Rr*Xs, Rr*Xe, l, sigmas2, sigmas_lines_2, sigmas2d_2, Sigmas2DLines, R_est_3*Rr', t_est_3);
            cnt = cnt+1;
        end
        for k = 1:size(R_est,3)
            R_est(:,:,k) = R_est(:,:,k)*Rr;
        end
        
    end
    
%     [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
    if mode == 1        
        return;
    end
    
    if size(R_est,3)>1
        [R_est_2, t_est_2] = dlsu(X, xn, Xs, Xe, l, Sigmas0, Sigmas3DLines, Sigmas2D_0, Sigmas2DLines, R_est(:,:,1), t_est(:,1));
        for k = 1:size(t_est_2,2)
            [R_est, t_est] = add_solution(R_est, t_est, R_est_2(:,:,k), t_est_2(:,k));
        end
    end
end