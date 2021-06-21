function [Rest, test] = upnp_pnpu(X, x, K, Sigmas3D, Sigmas2D, R_est_3, t_est_3)
    n = size(x,2);
    temp = K \ [x; ones(1,n)];
    I_norms = sqrt(sum(temp.*temp));
    I_normalized = temp ./ repmat(I_norms,3,1);    
    inds = 1:n;
    T_perturbed = [R_est_3 t_est_3];
     X = opengv('upnp',inds, X,I_normalized, T_perturbed);
    Rest = X(:,1:3);
    test = X(:,4);
end