function [XXw_s, XXw_e, Sigmas3DLines] = generate_3d_lines_with_covs(XXw_s, XXw_e, nln, nl_3d, isotropic)
    % shift 3D line endpoints
    shift_len = 0.1;
    dXw = XXw_e - XXw_s;                    
    nln_gen = size(XXw_s, 2);
    XXw_s = XXw_s + repmat(randn(1, nln_gen), 3, 1).*dXw*shift_len;
    XXw_e = XXw_e + repmat(randn(1, nln_gen), 3, 1).*dXw*shift_len;

    % add 3d noise
    nlsamples = ones(nln_gen, 1) / nln_gen;
    nnl = round(nln_gen * nlsamples);
    nls = zeros(3,nln_gen);
    id  = 1;
    %random variances x dim
    for idnl = 1:length(nlsamples)
        sigma = nl_3d(mod(idnl, length(nl_3d))+1);
        nls(1,id:id+nnl(idnl)-1) = sigma * ones(1,nnl(idnl));
        if(isotropic)
            nls(2,id:id+nnl(idnl)-1) = sigma * ones(1,nnl(idnl));
            nls(3,id:id+nnl(idnl)-1) = sigma * ones(1,nnl(idnl));
        else
            nls(2,id:id+nnl(idnl)-1) = sigma * rand(1,nnl(idnl));
            nls(3,id:id+nnl(idnl)-1) = sigma * rand(1,nnl(idnl));
        end
        id = id+nnl(idnl);
    end
    %noisy samples    
    Sigmas3DLines = zeros(6,6,nln);

    noise_3d_ln_inds = randperm(nln_gen);
    nls_lines = nls(:, noise_3d_ln_inds);

    for id = 1:length(nls(1,:))
        randR = rodrigues(randn(3,1));            
        sigmas = nls_lines(:, id);
        noise_start = randn(3, 1).*sigmas;
        randR_start = rodrigues(randn(3,1));
        Sigmas3DLines(1:3, 1:3, id) = randR_start * diag(sigmas.^2) * randR_start';
        XXw_s(:, id) = XXw_s(:, id) + randR_start * noise_start;
        noise_end = randn(3, 1).*sigmas;
        randR_end = rodrigues(randn(3,1));
        Sigmas3DLines(4:6, 4:6, id) = randR_end * diag(sigmas.^2) * randR_end';
        XXw_e(:, id) = XXw_e(:, id) + randR_end * noise_end;
    end
end