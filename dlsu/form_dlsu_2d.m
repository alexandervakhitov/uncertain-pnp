function A2 = form_dlsu_2d(X, x, Xs, Xe, l, Sigmas, SigmasLines, Sigmas2D, SigmasLines2D, R_est, t_est)
    np = size(X, 2);
    nl = size(Xs, 2);
    T = zeros(3, 3);
    Tp = zeros(2, 3, np);
    A = zeros(3, 9);    
    S2s = zeros(2, 2, np);
    Api = zeros(2, 9, np);
    
    depths_est = [];
    Xsc = [];
    Xec = [];
    if size(R_est, 1)>0
        Xc = R_est * X + repmat(t_est, 1, np);
        depths_est = Xc(3,:);    
        Xsc = R_est * Xs + repmat(t_est, 1, nl);
        Xec = R_est * Xe + repmat(t_est, 1, nl);    
    end
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
            %SigmaAdd = depths_est(i)*depths_est(i)*Sigmas2D(i)*Sigmas2D(i)*eye(2);%
            SigmaAdd = depths_est(i)*depths_est(i)*Sigmas2D(:,:,i);%
        end 
        SigmaUnc = SigmaAdd+1e-9*eye(2); %Sigma2 + 
        S2s(:,:,i) = inv(SigmaUnc);                
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
    Tl = zeros(3, 3, nl);
    Ali = zeros(2, 9, nl);
%     for i = 1:nl
%         Sigma2 = zeros(2,2);
%         Sigma2(1,1) = l(:, i)'*SigmasLines(1:3,1:3,i)*l(:, i);
%         Sigma2(1,2) = l(:, i)'*SigmasLines(1:3,4:6,i)*l(:, i);
%         Sigma2(2,2) = l(:, i)'*SigmasLines(4:6,4:6,i)*l(:, i);        
% %         if (size(R_est, 1) > 0)
% %             Sigma2(1,1) = Sigma2(1,1) + Xsc(:, i)'*SigmasLines2D(:,:,i)*Xsc(:, i);
% %             Sigma2(2,2) = Sigma2(2,2) + Xec(:, i)'*SigmasLines2D(:,:,i)*Xec(:, i);
% %             Sigma2(1,2) = Sigma2(1,2) + Xsc(:, i)'*SigmasLines2D(:,:,i)*Xec(:, i);
% %         end
%         Sigma2(2,1) = Sigma2(1,2);                        
%         S2ls(:, :, i) = inv(Sigma2);
%         
%         Ti = [l(:, i)'; l(:, i)'];
%         Xsi = Xs(:, i);
%         Asi = l(:, i)'*[Xsi', zeros(1,6);
%               zeros(1,3), Xsi', zeros(1,3);
%               zeros(1,6), Xsi'];
%         Xei = Xe(:, i);
%         Aei = l(:, i)'*[Xei', zeros(1,6);
%               zeros(1,3), Xei', zeros(1,3);
%               zeros(1,6), Xei'];
%         Ai = [Asi; Aei];
%         Ali(:, :, i) = Ai;
%         Tl(:, :, i) = Ti'*S2ls(:,:,i)*Ti;
%         T = T + Tl(:, :, i);                
%         A = A + Ti'*S2ls(:, :, i)*Ai;
%     end
    A1 = A;
    
    A2 = zeros(9,9);
    for i = 1:np                
        Ai = Api(:,:,i) - Tp(:, :, i)*(T\A1);
        Aim = Ai'*S2s(:, :, i)*Ai;        
        A2 = A2 + Aim;
    end
    for i = 1:nl        
        Ti = [l(:, i)'; l(:, i)'];
        Ai = Ali(:, :, i)-Ti*(T\A1);
        A2 = A2 + Ai'*S2ls(:, :, i)*Ai;        
    end
end