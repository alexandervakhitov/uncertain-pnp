function [R, t] = MLPNP_refine(R, t, points3D, v, cov)
    use_cov = 1;
    
    use_lsq = 0;
    
    nrPts = size(points3D,2);

    r = zeros(3,nrPts);
    s = zeros(3,nrPts);
    cov_reduced = zeros(2,2,nrPts);

    % test planarity, only works well if the scene is really planar
    % quasi-planar won't work very well
    S = points3D*points3D';
    
    planar = 0;
    
    if (rank(S) == 2)
        planar = 1;            
    end

    % compute null spaces of bearing vector v: null(v')
    for i=1:nrPts
        null_2d = null(v(1:3,i)');
        r(:,i) = null_2d(:,1);
        s(:,i) = null_2d(:,2);
        if use_cov
            tmp = reshape(cov(:,i),3,3);
            cov_reduced(:,:,i) = (null_2d'*tmp*null_2d)^-1;
        end
    end

    % stochastic model
    Kll = eye(2*nrPts,2*nrPts);
    % if (normalize)
    %     points3Dn = normc(points3Dn);
    % end    
    for i=1:nrPts
       if (use_cov)
            Kll(2*i-1:2*i,2*i-1:2*i) = cov_reduced(:,:,i);
       end
       if (use_lsq)
           Kll(2*i-1:2*i,2*i-1:2*i) = sqrtm(cov_reduced(:,:,i));
       end
       % r11           
    end



    if use_lsq        
        x = [Rodrigues2(R)', t']';
        fun = @(x) residualsAndJacobianLsq(x, r, s, points3D, Kll);
        OPTIONS = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'display', 'off', 'SpecifyObjectiveGradient', true);    %,'SpecifyObjectiveGradient',true
        x_fin = lsqnonlin(fun, x, zeros(0,0), zeros(0,0), OPTIONS);
        R = rodrigues(x_fin(1:3));
        t = x_fin(4:6);
    else
        optimFlags.epsP  = 1e-6;
        optimFlags.epsF  = 1e-6;
        optimFlags.maxit = 5;
        optimFlags.tau   = 1e-4;
        T = eye(4);
        T(1:3, 1:3) = R;
        T(1:3, 4) = t;
        [T, statistics] = optim_MLPnP_GN(T, points3D, r, s, Kll, optimFlags);
        R = T(1:3, 1:3);
        t = T(1:3, 4);
    end
end