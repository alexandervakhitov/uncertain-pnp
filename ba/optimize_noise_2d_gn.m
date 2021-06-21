function [R, t, statistic] = optimize_noise_2d_gn(R0, t0, XX, xx, sigmas_2d, optimFlags)
    [rvec, jac] = rodrigues(R0);
    x = [rvec; t0];    
    n = size(XX, 2);
    sigmas_2d_inv = zeros(2*n, 2*n);
    for i = 1:n        
        sigmas_2d_inv(2*i-1:2*i, 2*i-1:2*i) = inv(sigmas_2d(:, :, i));
    end    
    % redundancy
    redundanz = 2*n - length(x);    
    % optim params
    epsParam    = optimFlags.epsP;
    epsFunc     = optimFlags.epsF;

    % iteration params
    cnt = 0;
    stop = false;
    invKll = sigmas_2d_inv;
    while cnt < optimFlags.maxit && stop == 0
        [r, J] = residual_2d(x, XX, xx);                          
        % design matrix
        N = J.'*invKll*J;
        % System matrix
        g = J.'*invKll*r;

        dx = pinv(N)*g;
        if (max(abs(dx)) > 20 || min(abs(dx)) > 1)
            break;
        end
        dl = J*dx(1:end);

        if max(abs(dl)) < epsFunc || max(abs(dx(1:end))) < epsParam  
            x = x-dx;
            break;
        else
            % update parameter vector
            x = x-dx;              
        end  
        cnt = cnt+1;
    end % while loop


    % minimal to homogeneous
    R = rodrigues(x(1:3));
    t = x(4:6);

    % empirical variance factor
    resV = r.'*invKll*r;
    if redundanz > 0
        if redundanz < n
            s0 = 1;  
        else
            s0 = resV / redundanz;
        end
    else
        s0 = NaN;
    end
    % variance-covariance matrix
    Qxx = pinv(N);       
    % cofactor matrix of "adjusted observations"
    Qldld = J*Qxx*J';

    statistic = {resV, r, J, Qxx, s0, Qldld, sqrt(s0.*diag(Qxx))};

end

