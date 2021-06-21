function [R_est, t_est] = dlsuf_m_pnpu(X, x, K, Sigmas, Sigmas2D, R_est_3, t_est_3, mode, mode2)
    if nargin < 8
        mode = 1;
    end
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas0 = Sigmas;
    Sigmas2D_0 = Sigmas2D;
    if mode2 == 1
        [Sigmas, sigmas2] = replace_with_trace(Sigmas);
        [Sigmas2D, sigmas2d_2] = replace_with_trace(Sigmas2D);
    end
%     Xm = mean(X, 2);
    n = size(X, 2);    
    %Cp1 = mean(Xw, 2);
%     Cp1 = zeros(3,1);
%     S_sum = zeros(3,3);
%     for i = 1:n
%         Sci = inv(Sigmas(:,:,i)+1e-9*eye(3));
%         Cp1 = Cp1 + Sci*X(:, i);
%         S_sum = S_sum + Sci;
%     end
%     Si = inv(S_sum);
%     Xm = Si*Cp1;
    
%     for i = 1:n
%         Sigmas(:,:,i) = Sigmas(:,:,i)-Si;        
%     end
%     X = X - repmat(Xm, 1, n);
%     t_est_3 = t_est_3 + R_est_3*Xm;
    
%     xn = normalize_points(x, K);
%     Sigmas2D = Sigmas2D / K(1,1)/K(1,1);
    
    %[R_est, t_est] = dlsu(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3, t_est_3);
    xn = x;
    if mode2 == 1
        [R_est, t_est] = dlsu_fast(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), sigmas2, zeros(6,6,0), sigmas2d_2, zeros(3,3,0), R_est_3, t_est_3);
    else
        [R_est, t_est] = dlsu(X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3, t_est_3);
    end
    if (size(R_est, 3) == 0)
        cnt = 0;
        while cnt < 5 && size(R_est, 3) == 0
            Rr = rodrigues(3*randn(3,1));        
            %[R_est, t_est] = dlsu(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3*Rr', t_est_3);
            if mode2 == 1
                [R_est, t_est] = dlsu_fast(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), sigmas2, zeros(6,6,0), sigmas2d_2, zeros(3,3,0), R_est_3*Rr', t_est_3);
            else
                [R_est, t_est] = dlsu(Rr*X, xn, zeros(3,0), zeros(3,0), zeros(3, 0), Sigmas, zeros(6,6,0), Sigmas2D, zeros(3,3,0), R_est_3*Rr', t_est_3);
            end
            cnt = cnt+1;
        end
        for k = 1:size(R_est,3)
            R_est(:,:,k) = R_est(:,:,k)*Rr;
        end
        
    end
    
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end