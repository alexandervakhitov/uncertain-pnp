function [Rs, ts] = dlsl(X, x, Xs, Xe, l)
    np = size(X, 2);
    nl = size(Xs, 2);
    T = zeros(3, 3);
    Tp = zeros(2, 3, np);
    A = zeros(3, 9);    
    Api = zeros(2, 9, np);
    
    for i = 1:np    
        x1 = x(1,i);
        x2 = x(2,i);        
        Ti = [-1,0,x1; 
            0,-1,x2];        
        T = T + Ti'*Ti;
        Tp(:, :, i) = Ti;
        
        Xp = X(:, i);
        Ai = [-Xp', zeros(1,3), x1*Xp';
              zeros(1,3), -Xp', x2*Xp'];
        Api(:,:,i) = Ai;
        A = A + Ti'*Ai;
    end    
    nl = size(l, 2);
    Tl = zeros(3, 3, nl);
    Ali = zeros(2, 9, nl);
    for i = 1:nl               
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
        Tl(:, :, i) = Ti'*Ti;
        T = T + Tl(:, :, i);                
        A = A + Ti'*Ai;
    end
    A1 = A;
    
    A2 = zeros(9,9);
    for i = 1:np                
        Ai = Api(:,:,i) - Tp(:, :, i)*(T\A1);
        Aim = Ai'*Ai;        
        A2 = A2 + Aim;
    end
    for i = 1:nl        
        Ti = [l(:, i)'; l(:, i)'];
        Ai = Ali(:, :, i)-Ti*(T\A1);
        A2 = A2 + Ai'*Ai;        
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
    for i = 1:size(sols, 2)
        s_cand = sols(:, i);        
        if (isreal(s_cand))
            s = s_cand;
            Rc = (1-s'*s)*eye(3)+2*cpmat(s)+2*s*s';
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

