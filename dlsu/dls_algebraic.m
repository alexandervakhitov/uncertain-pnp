function [Rs, ts] = dls_algebraic(X, x, do_refine)
    if nargin < 3
        do_refine = 1;
    end
        
    
    tic;
    
    np = size(X, 2);        
    
    T = zeros(3, 3);
    Tp = zeros(2, 3, np);
    A = zeros(3, 9);    
    S2s = zeros(2, 2, np);
    Api = zeros(2, 9, np);
            
    for i = 1:np                
        x1 = x(1,i);
        x2 = x(2,i);
        
        S2s(:,:,i) = eye(2);
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
    A1 = A;
    
    A2 = zeros(9,9);
    for i = 1:np                
        Ai = Api(:,:,i) - Tp(:, :, i)*(T\A1);
        Aim = Ai'*S2s(:, :, i)*Ai;        
        A2 = A2 + Aim;
    end
    
    tcm = toc;
    
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
    
    tic;
    
    sols = solver_opt_pnp_hesch2(A3);
    
    tcs = toc;
    
    tic;
    
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
                
                sp = -s;
                if do_refine == 1
                    for it = 1:1                
                        g = dls_g(sp(1), sp(2), sp(3), A2);
                        H = dls_H(sp(1), sp(2), sp(3), A2);
                        sp = sp - H\g;                
                        Rit = 1.0/(1+sp'*sp)*(1-sp'*sp)*eye(3)+2*cpmat(sp)+2*sp*sp';  
                        r = Rit(:);
                        costval(it) = r'*A2*r;
                        sps(it,:) = sp;
                    end
                end
                sfin = -sp;
                Rc = 1.0/(1+sp'*sp)*((1-sfin'*sfin)*eye(3)+2*cpmat(sfin)+2*sfin*sfin');              
                Rct = Rc';
                r = Rct(:);

                tc = -T\A1*r;
                Rs(:,:,i) = Rc;
                ts(:,i) = tc;
            end                       
            
        end                
    end
   tcr = toc;
   Rs = Rs(:,:, valid_inds);
   ts = ts(:, valid_inds);
end