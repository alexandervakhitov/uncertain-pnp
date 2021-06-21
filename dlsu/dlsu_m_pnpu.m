function [R_est, t_est] = dlsu_m_pnpu(X, x, Sigmas, Sigmas2D, R_est_3, t_est_3, mode)
    if nargin < 7
        mode = 1;
    end
    
    tic;
    sigmas2 = compute_traces(Sigmas);
    sigmas2d_2 = compute_traces(Sigmas2D);        
    t_start = toc;
    
    tic;
    xn = x;       
    [R_est, t_est] = dlsu_fast(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), sigmas2, zeros(6,6,0), sigmas2d_2, zeros(3,3,0), R_est_3, t_est_3);
    if (size(R_est, 3) == 0)
        cnt = 0;
        while cnt < 5 && size(R_est, 3) == 0
            Rr = rodrigues(3*randn(3,1));        
            %[R_est, t_est] = dlsu(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3*Rr', t_est_3);
            if size(R_est_3, 1) > 0
                R_est_new = R_est_3*Rr';
            else
                R_est_new = R_est_3;            
            end
            t_est_new = t_est_3;
            [R_est, t_est] = dlsu_fast(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), sigmas2, zeros(6,6,0), sigmas2d_2, zeros(3,3,0), R_est_new, t_est_new);
            cnt = cnt+1;
        end
        for k = 1:size(R_est,3)
            R_est(:,:,k) = R_est(:,:,k)*Rr;
        end        
    end
    t_m1 = toc;
    
%     [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
    if mode == 1        
        t_fin = toc;
%         fprintf('dlsu %f %f | %f\n', t_start, t_m1, t_start+t_m1 );
        return;
    end
    
    tic;
    
    if size(R_est,3)>=1
        [R_est_2, t_est_2] = dlsu(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est(:,:,1), t_est(:,1));
        for k = 1:size(t_est_2,2)
            [R_est, t_est] = add_solution(R_est, t_est, R_est_2(:,:,k), t_est_2(:,k));
        end
    end
    
    t_fin = toc;
%     fprintf('dlsu2 %f %f %f | %f \n', t_start, t_m1, t_fin, t_fin+t_m1+t_start);
end