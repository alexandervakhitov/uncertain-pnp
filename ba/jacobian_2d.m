function J_p = jacobian_2d(XXc, xx_p, rot_jac, XX)
    n = size(XX, 2);    
    J_p = zeros(2, 6, n);       
%     J_p(1:2, 1:2, :) = repmat(reshape(eye(2), 2, 2, 1), 1, 1, n);
    J_p(1, 1, :) = reshape(1./XXc(3, :), 1, 1, n);
    J_p(2, 2, :) = reshape(1./XXc(3, :), 1, 1, n);
    J_p(1:2, 3, :) = reshape(-xx_p ./ repmat(XXc(3, :), 2, 1), 2, 1, n);
    rot_jac = reshape(rot_jac, 3, 3, 3);
    J_3d = zeros(3, 6, n);
    J_3d(1, 4, :) = 1;
    J_3d(2, 5, :) = 1;
    J_3d(3, 6, :) = 1;
    for rpi = 1:3
       J_3d(:, rpi, :) = reshape(rot_jac(:, :, rpi) * XX, 3, 1, n);
    end
    for i = 1:n           
       J_p(:, :, i) = J_p(:, 1:3, i) * J_3d(:, :, i);           
    end
end