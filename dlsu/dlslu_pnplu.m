function [R_est, t_est] = dlslu_pnplu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3, Xs, Xe, l, Sigmas3DLines, Sigmas2DLines, xs, xe)
%poses = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    Sigmas = replace_with_trace(Sigmas3D);
    SigmasLines = Sigmas3DLines;
    SigmasLines(1:3,1:3,:) = replace_with_trace(Sigmas3DLines(1:3,1:3,:));
    SigmasLines(1:3,4:6,:) = replace_with_trace(Sigmas3DLines(1:3,4:6,:));
    SigmasLines(4:6,4:6,:) = replace_with_trace(Sigmas3DLines(4:6,4:6,:));
    SigmasLines(4:6,1:3,:) = replace_with_trace(Sigmas3DLines(4:6,1:3,:));
% 
    
%     xn = normalize_points(x, K);
%     ln = normalize_lines(l, K);   
%     X = zeros(3,0);
%     xn = zeros(2,0);
%     Sigmas2D  = zeros(2,2,0);
%     Sigmas = zeros(3,3,0);
%     dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)

    Xd = Xe - Xs;
    if (size(Xd, 2) > 0)
        Xdn = sqrt(Xd(1, :).^2 + Xd(2, :).^2 + Xd(3, :).^2);
        Xd = Xd ./ repmat(Xdn, 3, 1);
        Xe = Xs + Xd;
    end

    [R_est, t_est] = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, Sigmas2DLines, R_est_3, t_est_3);
    if (size(R_est, 3) == 0)
        cnt = 0;
        while cnt < 5 && size(R_est, 3) == 0
            Rr = rodrigues(3*randn(3,1));        
            [R_est, t_est] = dlsu(Rr*X, x, Rr*Xs, Rr*Xe, l, Sigmas, SigmasLines, Sigmas2D, Sigmas2DLines, R_est_3*Rr', t_est_3);
            cnt = cnt+1;
        end
        for k = 1:size(R_est,3)
            R_est(:,:,k) = R_est(:,:,k)*Rr;
        end
        
    end
    [R_est, t_est] = add_solution(R_est, t_est, R_est_3, t_est_3);
    
end
