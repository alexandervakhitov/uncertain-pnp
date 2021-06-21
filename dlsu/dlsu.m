function [Rs, ts, costs] = dlsu(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est, do_refine)

    if nargin < 12
        do_refine = 1;
    end                   
    
    np = size(X, 2);        
    
    nl = size(Xs, 2);
    T = zeros(3, 3);
    Tp = zeros(2, 3, np);
    A = zeros(3, 9);    
    S2s = zeros(2, 2, np);
    Api = zeros(2, 9, np);
    
    
    Xsc = [];
    Xec = [];
    if size(R_est, 1)>0
        Xc = R_est * X + repmat(t_est, 1, np);
        depths_est = Xc(3,:);    
        Xsc = R_est * Xs + repmat(t_est, 1, nl);
        Xec = R_est * Xe + repmat(t_est, 1, nl);    
        depths_start_est = Xsc(3, :);
        depths_end_est = Xec(3, :);
    else
        depths_est = ones(1, np);
        depths_start_est = ones(1, nl);
        depths_end_est = ones(1, nl);
        R_est = eye(3);
    end
    
    depths_rep = repmat(reshape(depths_est, 1, 1, np), 2, 2, 1);
    SigmaAddAll = depths_rep .* depths_rep .* Sigmas2D;
    for i = 1:np        
        Sigma2 = zeros(2,2);
        x1 = x(1,i);
        x2 = x(2,i);
        Sigma3 = R_est * Sigmas(:,:,i) * R_est';
        Sigma2(1,1) = Sigma3(3,3)*x1*x1-2*Sigma3(1,3)*x1+Sigma3(1,1);
        Sigma2(2,2) = Sigma3(3,3)*x2*x2-2*Sigma3(2,3)*x2+Sigma3(2,2);
        Sigma2(1,2) = Sigma3(3,3)*x1*x2-Sigma3(1,3)*x2-Sigma3(2,3)*x1+Sigma3(1,2);
        Sigma2(2,1) = Sigma2(1,2);
        SigmaAdd = zeros(2,2);
        if size(R_est, 1)>0
            SigmaAdd = SigmaAddAll(:, :, i);
        end 
        SigmaUnc = Sigma2 + SigmaAdd;
        
        S2s(:,:,i) = inv(SigmaUnc + 1e-6 * eye(2)); 
        Ti = [-1,0,x1; 
            0,-1,x2];        
        T = T + Ti'*S2s(:,:,i)*Ti;
        Tp(:, :, i) = Ti;
        
        Xp = X(:, i);
        Ai = [-Xp', zeros(1,3), x1*Xp';
              zeros(1,3), -Xp', x2*Xp'];
        Api(:,:,i) = Ai;
        A = A + Ti'*S2s(:, :, i)*Ai;
    end
    S2ls = zeros(2, 2, nl);
    nl = size(l, 2);
    Tl = zeros(2, 3, nl);
    Ali = zeros(2, 9, nl);
    for i = 1:nl
        Sigma2 = zeros(2,2);
        Sigma2(1,1) = l(:, i)'*R_est*SigmasLines(1:3,1:3,i)*R_est'*l(:, i);
        Sigma2(1,2) = l(:, i)'*R_est*SigmasLines(1:3,4:6,i)*R_est'*l(:, i);
        Sigma2(2,2) = l(:, i)'*R_est*SigmasLines(4:6,4:6,i)*R_est'*l(:, i);        
        Sigma2(2,1) = l(:, i)'*R_est*SigmasLines(4:6,1:3,i)*R_est'*l(:, i);
        if (size(R_est, 1) > 0)
            Sigma2 = Sigma2 + diag([depths_start_est(1, i).^2, depths_end_est(1, i).^2]) * SigmasLines2D(i);
        end
        S2ls(:, :, i) = inv(Sigma2 + 1e-9*eye(2));
        
        Ti = [l(:, i)'; l(:, i)'];
        Xsi = Xs(:, i);
        Asi = l(:, i)'*[Xsi', zeros(1,6);
              zeros(1,3), Xsi', zeros(1,3);
              zeros(1,6), Xsi'];
        Xei = Xe(:, i);
        Aei = l(:, i)'*[Xei', zeros(1,6);
              zeros(1,3), Xei', zeros(1,3);
              zeros(1,6), Xei'];
        Ai = [Asi; Aei];
        Ali(:, :, i) = Ai;        
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
        s_cand = sols(:, i);        
        if (isreal(s_cand))
            s = s_cand;
            Rc = rot_from_sol(s);
            Rct = Rc';
            r = Rct(:);                                    
            tc = -T\A1*r;
            Rs(:, :, i) = Rc;
            ts(:, i) = tc;
            % check that the points are infront of the center of perspectivity            
            Xc = Rc*X + tc;
                        
            if isempty(find(Xc(3,:) < 0))
                
                sp = -s;
                sinit = sp;
                if do_refine == 1
                    Rit = rot_from_sol(sp);  
                    r = Rit(:);
                    costval(1) = r'*A2*r;
                    for it = 1:2                
                        g = dls_g(sp(1), sp(2), sp(3), A2);
                        H = dls_H(sp(1), sp(2), sp(3), A2);                        
                        sp = sp - H\g;                                        
                        Rit = rot_from_sol(sp);  
                        r = Rit(:);
                        costval(it+1) = r'*A2*r;                        
                        sps(it,:) = sp;
                    end                    
                    if (costval(3) > costval(1))
                        sp = sinit;
                    end
                end                
                sfin = -sp;
                Rc = rot_from_sol(sfin);              
                Rct = Rc';
                r = Rct(:);                                
                tc = -T\A1*r;
                Rs(:,:,i) = Rc;
                ts(:,i) = tc;
                
                Xce = Rc * X + repmat(tc, 1, np);
                xce = Xce(1:2, :) ./ repmat(Xce(3, :), 2, 1);
                fin_cost = max(sqrt(sum((xce - x).^2, 1)));
                valid_inds = [valid_inds; i];
                costs = [costs; fin_cost];
            end                                   
        end                
    end
   tcr = toc;
   Rs = Rs(:,:, valid_inds);
   ts = ts(:, valid_inds);
end

function R = rot_from_sol(sp)
    R = 1.0/(1+sp'*sp)*((1-sp'*sp)*eye(3)+2*cpmat(sp)+2*sp*sp');  
end