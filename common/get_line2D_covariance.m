function Sigma_theory = get_line2D_covariance(Sigma_xsu, Sigma_xeu, xs0h, xe0h)
    l0_h = cross(xs0h, xe0h);
    l0_h_norm = norm(l0_h(1:2));
    Sigma_xsu_h = zeros(3,3);
    Sigma_xsu_h(1:2, 1:2) = Sigma_xsu;    
    Sigma_xeu_h = zeros(3,3);
    Sigma_xeu_h(1:2, 1:2) = Sigma_xeu;    
    Sigma_l_h = ( cpmat(xe0h) * Sigma_xsu_h * cpmat(xe0h)' + ...
        cpmat(xs0h) * Sigma_xeu_h * cpmat(xs0h)');

    J = zeros(3,3);
    l0_h_2 = l0_h(1:2);
    J(1:2,1:2) = 1.0/l0_h_norm * eye(2) - 1.0/l0_h_norm^3*(l0_h_2*l0_h_2');
    J(3, 1:2) = -1.0/l0_h_norm^3*l0_h(3)*l0_h_2';
    J(3, 3) = 1/l0_h_norm;
    
%     J = eye(3);
    Sigma_theory = J * Sigma_l_h *J';
end