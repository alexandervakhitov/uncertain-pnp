function [R_est, t_est] = robust_dls_algebraic(X, x)
    [R_est, t_est] = dls_algebraic(X, x);
    if (size(R_est, 3) == 0)
        cnt = 0;
        while cnt < 5 && size(R_est, 3) == 0
            randax = 3*randn(3,1);
            randax(randi(3)) = 0;
            Rr = rodrigues(randax);            
            [R_est, t_est] = dls_algebraic(Rr*X, x);            
            for k = 1:size(R_est,3)
                R_est(:,:,k) = R_est(:,:,k)*Rr;
            end                        
            cnt = cnt+1;
        end
    end
end