function [r, J] = residual_2d(x, XX, xx, Xs, Xe, ll, sigmas_2d_inv_sqrt, sigmas_2d_lines_inv_sqrt)
%     if nargin < 4
%         Xs = zeros(3, 0);
%         Xe = zeros(3, 0);
%         ll = zeros(3, 0);
%         sigmas_2d_lines_inv_sqrt = [];
%     end
    [R, rot_jac] = rodrigues(x(1:3));
    t = x(4:6);
    n = size(XX, 2);
    XXc = R * XX + repmat(t, 1, n);
    xx_p = XXc(1:2, :) ./ repmat(XXc(3, :), 2, 1);
    dx_w = zeros(size(xx));
    dx = xx_p - xx;
    if nargin > 6
        for i = 1:n
            dx_w(:, i) = sigmas_2d_inv_sqrt(:, :, i) * dx(:, i);
        end
    else        
        dx_w = dx;
    end
        
    nl = size(Xs, 2);
    Xsc = R * Xs + repmat(t, 1, nl);
    xxs_p = Xsc(1:2, :) ./ repmat(Xsc(3, :), 2, 1);
    dxs = sum(ll(1:2, :).*xxs_p, 1) + ll(3, :);
    Xec = R * Xe + repmat(t, 1, nl);
    xxe_p = Xec(1:2, :) ./ repmat(Xec(3, :), 2, 1);
    dxe = sum(ll(1:2, :).*xxe_p, 1) + ll(3, :);    
    if nargin > 6
        l_weight = reshape(sigmas_2d_lines_inv_sqrt, 1, nl);
    else
        l_weight = ones(1, nl);
    end
    dl_w = [dxs .* l_weight; dxe .* l_weight ];        

    r = [dx_w(:); dl_w(:)];
    if nargout > 1   % Two output arguments
       nf = n + nl;
       J = zeros(2*nf, 6);
       J_p = jacobian_2d(XXc, xx_p, rot_jac, XX);
       if nargin > 6
           for i = 1:n
               J(2*i-1:2*i, :) = sigmas_2d_inv_sqrt(:, :, i) * J_p(:, :, i);
           end
       else
           for i = 1:n
               J(2*i-1:2*i, :) = J_p(:, :, i);
           end
       end
       l_rep = repmat(reshape(ll(1:2, :), 2, 1, nl), 1, 6, 1);
       J_s = jacobian_2d(Xsc, xxs_p, rot_jac, Xs);
       J_s = reshape(sum(J_s .* l_rep, 1), 6, nl)';
       J_e = jacobian_2d(Xec, xxe_p, rot_jac, Xe);
       J_e = reshape(sum(J_e .* l_rep, 1), 6, nl)';
       if nargin > 6
           for i = 1:nl
               J(2*n + 2*i - 1, :) = sigmas_2d_lines_inv_sqrt(i) * J_s(i, :);
               J(2*n + 2*i, :) = sigmas_2d_lines_inv_sqrt(i) * J_e(i, :);
           end
       else
           for i = 1:nl
               J(2*n + 2*i - 1, :) = J_s(i, :);
               J(2*n + 2*i, :) = J_e(i, :);
           end
       end
    end
end