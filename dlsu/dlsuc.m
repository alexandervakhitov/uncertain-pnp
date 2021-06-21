function [Rs, ts] = dlsuc(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    A2_2d = form_dlsu_2d(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est);
    A2_3d = form_dlsu_3d(X, x, Xs, Xe, l, Sigmas, SL, Sigmas2D, R_est_3, t_est_3);
    %QR : [1 s1 s2 s3 s1^2 s1s2 s1s3 s2^2 s2s3 s3^2] -> r
    QR = zeros(9, 10);
    QR(1, :) = [1,0,0,0, 1,0,0,-1, 0, -1];
    QR(2, :) = [0,0,0,-2, 0,2,0,0, 0, 0];
    QR(3, :) = [0,0,2,0, 0,0,2,0, 0, 0];
    QR(4, :) = [0,0,0,2, 0,2,0,0, 0, 0];
    QR(5, :) = [1,0,0,0, -1,0,0,1, 0, -1];
    QR(6, :) = [0,-2,0,0, 0,0,0,0, 2, 0];
    QR(7, :) = [0,0,-2,0, 0,0,2,0, 0, 0];
    QR(8, :) = [0,2,0,0, 0,0,0,0, 2, 0];
    QR(9, :) = [1,0,0,0, -1,0,0,-1, 0, 1];
    
    A3 = QR'*(A2_2d+A2_3d)*QR;
    
    sols = solver_opt_pnp_hesch2(A3);
    
    Rs = zeros(3, 3, length(sols));
    ts = zeros(3, length(sols));
    valid_inds = [];
    for i = 1:size(sols, 2)
        s_cand = sols(:, i);        
        if (isreal(s_cand))
            s = s_cand;
            Rc = (1-s'*s)*eye(3)+2*pnp3d.cpmat(s)+2*s*s';
            Rc = Rc / (det(Rc))^(1/3);
            Rct = Rc';
            r = Rct(:);
            tc = -T\A1*r;
            Rs(:,:,i) = Rc;
            ts(:,i) = tc;
            % check that the points are infront of the center of perspectivity            
            Xc = Rc*X + tc;
            if isempty(find(Xc(3,:) < 0))
                valid_inds = [valid_inds; i];
            end           
            
        end                
    end
   Rs = Rs(:,:, valid_inds);
   ts = ts(:, valid_inds);
end

