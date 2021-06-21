function [XXw, Sigmas3D] = generate_3d_points_with_covs(XXw, nl, nlsamples, isotropic )
    npt = size(XXw, 2);
    nnl = round(npt * nlsamples);
    nls = zeros(3,npt);
    id  = 1;
    %random variances x dim
    for idnl = 1:length(nl)
        sigma = nl(idnl);
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
    noise = zeros(3,1);
    Sigmas3D = zeros(3, 3, npt);
    for id = 1:length(nls(1,:))
        sigmas = nls(:,id);
        noise(1) = randn * sigmas(1);
        noise(2) = randn * sigmas(2);
        noise(3) = randn * sigmas(3);

        %randR = eye(3);
        randR = rodrigues(randn(3,1));
        Sigmas3D(:,:,id) = randR * diag(sigmas.^2) * randR';            
        XXw(:,id) = XXw(:,id) + randR * noise;
    end
end