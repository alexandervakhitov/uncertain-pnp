function [Rest, test] = mlpnp_pnpu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3)
    n = size(X, 2);
    v = zeros(3, n);
    SigmasBear = zeros(9, n);
    for i = 1:n
        xh = K\[x(:, i); 1];
        v(:, i) = xh/norm(xh);
        S2d = Sigmas2D(:,:,i);
        S2dh = zeros(3,3);
        S2dh(1:2,1:2) = S2d;
        J = 1/norm(xh(1:2))*(eye(3) - v(:, i)*v(:, i)');
        Sbear = J*S2dh*J';
        SigmasBear(:, i) = reshape(Sbear, 9, 1);
    end
    [T, stats] = MLPnP(X, v, SigmasBear);
    Rest = T(1:3, 1:3);
    test = T(1:3, 4);
end