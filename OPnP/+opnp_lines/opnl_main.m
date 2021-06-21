function [R0 t0 error0 flag as] = opnl_main(Xs, Xe, l2d, as)
    if (nargin < 4)
        as = 0;            
    end
    if (as == 2)
        Rfix = opnp_lines.rodrigues([1 0 0]');
    end
    if (as == 1)
%         Rfix = rodrigues(randn(3,1));
        Rfix = opnp_lines.rodrigues([1 1 1]');
    end
    
    Xso = Xs;
    Xeo = Xe;    

    if (as > 0)
        Xs = Rfix*Xs;
        Xe = Rfix*Xe;
    end
    
    L = l2d';
    Xs1 = Xs';
    Xe1 = Xe';
    Lrep = repmat(L, 1, 3);
    N011 = Lrep.*[repmat(Xs1(:, 1), 1, 3) repmat(Xs1(:, 2), 1, 3) repmat(Xs1(:, 3), 1, 3)];
    N012 = Lrep.*[repmat(Xe1(:, 1), 1, 3) repmat(Xe1(:, 2), 1, 3) repmat(Xe1(:, 3), 1, 3)];
    
    N0 = [N011; N012];
    RtoPol  = rToPol();
    RtoPol = RtoPol(2:end, :);    

    Lb = repmat(L, 2, 1);
    
%     [U,S,V] = svd(Lb);
%     U1 = U(:, 1:3)';
%     LbTN0 = V*S(1:3,1:3)*U1*N0;
%     Li = inv(V*S(1:3,1:3).^2*V');
%     M0 = N0'*N0 - LbTN0'*Li*LbTN0;
    M0 = N0'*N0 - N0'*Lb*inv(Lb'*Lb)*Lb'*N0;
%     M0 = N0'*N0 - N0'*(U(:, 1:3)*U( :, 1:3)')*N0;
%     LiLbtN0 =  V*inv(S(1:3,1:3))*U(:, 1:3)'*N0;
    LiLbtN0 = inv(Lb'*Lb)*Lb'*N0;
%     Pm = Lb*P0;    
%     M0 = N0'*N0 - N0'*Pm*N0;
    M = RtoPol*M0*RtoPol';
    t0 = tic;
    [be,ce,de] = opnp_lines.solver_groebner_test(M(:));
    t1 = toc(t0);
    t0 = tic;
    solNum = length(be);
    Rs = zeros(3,3,solNum);
    ts = zeros(3, solNum);
    errs = zeros(solNum, 1);
%     Q = [A L]'*[A L];
    for i = 1:length(be)        
        b = be(i);
        c = ce(i);
        d = de(i);
        x = [1, b, c, d, b^2, b*c, b*d, c^2, c*d, d^2];        
        r = x*RtoPol;                
        Rr0 = reshape(r, 3, 3);
        scaleFactor = 1/ det(Rr0)^(1/3);
        Rr0 = Rr0*scaleFactor;   
        
        r0 = Rr0(:);
        fmin = r0'*M0*r0;
        
        [xopt, fminf] = opnp_lines.refine_pnpl(Rr0, M0, 10);
%         [xopt, fminf] = opnp_lines.refine_pnpl_levmar(Rr0, M0);
        
        if (fminf < fmin)
            Rr = opnp_lines.rodrigues(xopt);
            r = Rr(:);   
            fmin = fminf;
        else
            r = Rr0(:);
            Rr = Rr0;
        end
%         tr = -LiLbtN0*r;
%         if (as >0)
%             Rr = Rr*Rfix;
%         end
%         
%         if (tr(3) < 0)
%            Rr = [-Rr(:,1) -Rr(:, 2) Rr(:, 3)];
%            tr = -tr;
%         end
        Rs(:, :, i) = Rr;
%         ts(:, i) = tr;                
        errs(i) = fmin;%check_reproj(Rr, tr, l2d, Xso, Xeo);
    end
    t2 = toc(t0);
    error0 = min(abs(errs));
    del = 1e-4;
    rr = find(errs <= error0+del); 
    R0 = Rs(:,:,rr);
    t0 = ts(:,rr);    
    errs = errs(rr);
    
    for i = 1:length(rr)
        Rr0 = R0(:, :, i);
        r0 = Rr0(:);
        fmin = r0'*M0*r0;
        [xopt, fminf] = opnp_lines.refine_pnpl(Rr0, M0, 10);
%         [xopt, fmin] = refine_pnpl_levmar(Rr0, M0);
        
        if (fminf < fmin)
            Rr = opnp_lines.rodrigues(xopt);
            r = Rr(:);   
            fmin = fminf;
        else
            r = Rr0(:);
            Rr = Rr0;
        end
        tr = -LiLbtN0*r;
        if (as >0)
            Rr = Rr*Rfix;
        end
        
        if (tr(3) < 0)
           Rr = [-Rr(:,1) -Rr(:, 2) Rr(:, 3)];
           tr = -tr;
        end        
        
        R0(:, :, i) = Rr;
        t0(:, i) = tr;          
    end
    
    flag = 0;
    as = zeros(length(rr), 1);
    for i = 1:length(rr)
        if (length(rr) > 1)
            q = opnp_lines.rot2quat(R0(:, :, i));
        else
            q = opnp_lines.rot2quat(R0);
        end;
        if (abs(q(1)) < 0.01)% || errs(i) > 0.05)
            as(i) = 1;
        end
    end
%     fprintf('time %f %f\n', t1, t2);
end


function RtoPol  = rToPol
    RtoPol = zeros(11, 9);
    RtoPol(2, 1) = 1;
    RtoPol(6, 1) = 1;
    RtoPol(9, 1) = -1;
    RtoPol(11, 1) = -1;

    RtoPol(7, 2) = 2;
    RtoPol(5, 2) = -2;

    RtoPol(4, 3) = 2;
    RtoPol(8, 3) = 2;

    RtoPol(5, 4) = 2;
    RtoPol(7, 4) = 2;

    RtoPol(2, 5) = 1;
    RtoPol(6, 5) = -1;
    RtoPol(9, 5) = 1;
    RtoPol(11, 5) = -1;

    RtoPol(3, 6) = -2;
    RtoPol(10, 6) = 2;

    RtoPol(4, 7) = -2;
    RtoPol(8, 7) = 2;

    RtoPol(3, 8) = 2;
    RtoPol(10, 8) = 2;

    RtoPol(2, 9) = 1;
    RtoPol(6, 9) = -1;
    RtoPol(9, 9) = -1;
    RtoPol(11, 9) = 1;
    RtoPol(1,:) = -1;
end
function rep = check_reproj(R, t, l2d, Xs, Xe)
    tcop = repmat(t, 1, size(Xs, 2));
    Xs_cam = R*Xs + tcop;
    Xe_cam = R*Xe + tcop;
    Xs_cam = Xs_cam ./ repmat(Xs_cam(3, :), 3, 1);
    Xe_cam = Xe_cam ./ repmat(Xe_cam(3, :), 3, 1);
    reps1 = mean(abs(l2d(1,:).*Xs_cam(1,:) + l2d(2,:).*Xs_cam(2,:) + l2d(3,:).*Xs_cam(3,:)));
    reps2 = mean(abs(l2d(1,:).*Xe_cam(1,:) + l2d(2,:).*Xe_cam(2,:) + l2d(3,:).*Xe_cam(3,:)));    
    rep = 0.5*(reps1 + reps2);
end