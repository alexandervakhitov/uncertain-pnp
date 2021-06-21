function [Rs, ts, costs] = dlsu_fast(X, x, Xs, Xe, l, sigmas2, SigmasLines, ...
                                     sigmas2D_2, SigmasLines2D, R_est, t_est)
    if size(l, 1) == 0
        l = zeros(3, 0);
        Xs = zeros(3, 0);
        Xe = zeros(3, 0);        
    end
        
    np = size(X, 2);            
    
    nl = size(Xs, 2);
    T = zeros(3, 3);    
    A = zeros(3, 9);    
    S2s = zeros(2, 2, np);
    Api = zeros(2, 9, np);
    Tp = zeros(2, 3, np);
        
    if size(R_est, 1)>0
        Xc = R_est * X + repmat(t_est, 1, np);
        depths_est = Xc(3,:);                    
        Xsc = R_est * Xs + repmat(t_est, 1, nl);
        Xec = R_est * Xe + repmat(t_est, 1, nl);    
        start_depths_est = Xsc(3, :);
        end_depths_est = Xec(3, :);        
    else
        depths_est = ones(1, np) * t_est;        
        start_depths_est = ones(1, nl) * t_est;
        end_depths_est = ones(1, nl) * t_est;
    end
    sigmas_2d_add = depths_est.*depths_est.*sigmas2D_2;
    
    Tp(1:2, 1:2, :) = -repmat(reshape(eye(2), 2, 2, 1), 1, 1, np);
    Tp(1:2, 3, :) = reshape(x, 2, 1, np);
    Api(1, 1:3, :) = -reshape(X, 1, 3, np);
    Api(1, 7:9, :) = reshape(X.*repmat(x(1, :), 3, 1), 1, 3, np);    
    Api(2, 4:6, :) = -reshape(X, 1, 3, np);
    Api(2, 7:9, :) = reshape(X.*repmat(x(2, :), 3, 1), 1, 3, np);
    
    for i = 1:np        
        SigmaUnc = (sigmas_2d_add(i) + sigmas2(i)) * eye(2) + ...
                   sigmas2(i) * x(:, i) * x(:, i)';
        SigmaUncInv = inv(SigmaUnc + 1e-6 * eye(2));
        S2s(:,:,i) = SigmaUncInv;
        Ti = Tp(:, :, i);
        T = T + Ti'*S2s(:,:,i)*Ti;
        Ai = Api(:, :, i);
        A = A + Ti'*S2s(:, :, i)*Ai;
    end
    S2ls = zeros(2, 2, nl);    
    Tl = repmat(reshape(l, 1, 3, nl), 2, 1, 1);
    Ali = zeros(2, 9, nl);
    Xs1 = repmat(l(1, :), 3, 1) .* Xs;
    Xs2 = repmat(l(2, :), 3, 1) .* Xs;
    Xs3 = repmat(l(3, :), 3, 1) .* Xs;
    Ali(1, :, :) = reshape([Xs1; Xs2; Xs3], 1, 9, nl);
    Xe1 = repmat(l(1, :), 3, 1) .* Xe;
    Xe2 = repmat(l(2, :), 3, 1) .* Xe;
    Xe3 = repmat(l(3, :), 3, 1) .* Xe;
    Ali(2, :, :) = reshape([Xe1; Xe2; Xe3], 1, 9, nl);
    
    Sigmas_2D_Lines_start = start_depths_est.^2.*SigmasLines2D;
    Sigmas_2D_Lines_end = end_depths_est.^2.*SigmasLines2D;
    l_norms_sq = sum(l.^2, 1);
    for i = 1:nl
        Sigma2 = eye(2) * l_norms_sq(i) * SigmasLines(i);
        Sigma2 = Sigma2 + diag([Sigmas_2D_Lines_start(i), Sigmas_2D_Lines_end(i)]);

        S2ls(:, :, i) = inv(Sigma2 + 1e-9*eye(2));
        Ti = Tl(:, :, i);
        Ai = Ali(:, :, i);
        T = T + Ti'*S2ls(:,:,i)*Ti;                
        A = A + Ti'*S2ls(:, :, i)*Ai;
        Tl(:, :, i) = Ti;
    end
    A1 = A;
    
    A2 = zeros(9,9);
    for i = 1:np                
        Ai = Api(:,:,i) - Tp(:, :, i)*(T\A1);
        Aim = Ai'*S2s(:, :, i)*Ai;        
        A2 = A2 + Aim;
    end
    for i = 1:nl                
        Ai = Ali(:, :, i)-Tl(:,:,i)*(T\A1);
        A2 = A2 + Ai'*S2ls(:, :, i)*Ai;        
    end
        
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
    
    A3 = QR'*A2*QR;
        
    sols = solver_opt_pnp_hesch2(A3);
    
    Rs = zeros(3, 3, length(sols));
    ts = zeros(3, length(sols));
    valid_inds = [];
    costs = [];
    for i = 1:size(sols, 2)
        s = sols(:, i);        
        if (isreal(s))            
            Rc = 1.0/(1+s'*s)*((1-s'*s)*eye(3)+2*cpmat(s)+2*s*s');              
            Rct = Rc';
            r = Rct(:);
                                    
            tc = -T\A1*r;
            Rs(:,:,i) = Rc;
            ts(:,i) = tc;
            % check that the points are in front of the camera
            Xc = Rc*X + tc;
            if isempty(find(Xc(3,:) < 0))                                
                xc = Xc(1:2, :) ./ repmat(Xc(3, :), 2, 1);
                fin_cost = max(sqrt(sum((xc - x).^2, 1)));
                valid_inds = [valid_inds; i];
                costs = [costs; fin_cost];
            end                                   
        end                
    end
    
   Rs = Rs(:,:, valid_inds);
   ts = ts(:, valid_inds);
end