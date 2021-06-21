function A2 = form_dlsu_3d(X, x, Xs, Xe, l, Sigmas, SL, Sigmas2D, R_est_3, t_est_3)
    np = size(X, 2);
  
    T = zeros(3, 3);
    Tp = zeros(3, 3, np);
    A = zeros(3, 9);    
    for i = 1:np
        xp = [x(:, i); 1];
        J = 1/norm(xp)*eye(3)-1/norm(xp)^3*xp*xp';
        xp = xp/norm(xp);        
        Xp = X(:, i);
        Xpc = R_est_3*Xp + t_est_3;
        depth = Xpc(2);
        %S = R_est_3*inv(Sigmas(:,:,i)+1e-5*eye(3))*R_est_3';
        Sigma3d  = Sigmas(:,:,i) + 1e-5*eye(3);
        
%         Sigma3d(1:2,1:2) = Sigma3d(1:2,1:2) + depth*depth*Sigmas2D(:,:,i);
        Su = zeros(3);
        Su(1:2,1:2) = Sigmas2D(:,:,i);
        Su = J*Su*J';
        U = xp*xp';
        XM = Xpc*Xpc';        
        Sigma3d_from_2d = Su*(xp'*XM*xp) + U*XM*Su + Su*XM*U+(Xpc'*Su*Xpc)*U + Su*XM+XM*Su;
        Sigma3d = Sigma3d;% + Sigma3d_from_2d;
        S = R_est_3*inv(Sigma3d)*R_est_3';
        Ui = eye(3) - 1/(xp'*S*xp)*(xp*xp')*S;
        USU = Ui'*S*Ui;        
        Tp(:, :, i) = USU;
        T = T + USU;
        
        Ai = [Xp', zeros(1,6);
              zeros(1,3), Xp', zeros(1,3);
              zeros(1,6), Xp'];
        A = A + USU*Ai;            
    end
    nl = 0;%size(l, 2);
    Tl = zeros(3, 3, nl, 2);
    for i = 1:nl
        li = l(:, i) / norm(l(:, i));
        for j = 1:2
            Li = li*li';
            LSL = Li*SL(:, :, i, j)*Li;
            T = T + LSL;
            Tl(:, :, i, j) = LSL;
            Xl = zeros(3,1);
            if (j == 1)
                Xl = Xs(:, i);
            else
                Xl = Xe(:, i);
            end
            Ai = [Xl', zeros(1,6);
              zeros(1,3), Xl', zeros(1,3);
              zeros(1,6), Xl'];
            A = A + LSL*Ai;
        end
    end
    A1 = A;
    
    A2 = zeros(9,9);
    for i = 1:np
        Xp = X(:, i);
        Ai = [Xp', zeros(1,6);
              zeros(1,3), Xp', zeros(1,3);
              zeros(1,6), Xp'];
        Aim = Ai - T\A1;
        A2 = A2 + Aim'*Tp(:, :, i)*Aim;
    end
    for i = 1:nl
        for j = 1:2
            if (j == 1)
                Xl = Xs(:, i);
            else
                Xl = Xe(:, i);
            end
            Ai = [Xl', zeros(1,6);
              zeros(1,3), Xl', zeros(1,3);
              zeros(1,6), Xl'];
            Aim = Ai - T\A1;
            A2 = A2 + Aim'*Tl(:, :, i, j)*Aim;
        end
    end
end