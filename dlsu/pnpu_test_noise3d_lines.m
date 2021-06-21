function pnpu_test_noise3d_lines
    addpath('../');
    addpath('../dls_pnp_matlab');
    addpath('../p3p_code_final');
    addpath('../EPnP_matlab/EPnP');
    addpath('/home/alexander/materials/pnp3d/pnp3d/18');

    method_list = {@p3p_pnplu, @dls_pnplu, @dlsu_pnplu, @dlsl_pnplu, @dlslu_pnplu};   
%     method_list = {@p3p_pnplu, @dls_pnplu, @dlsu_pnplu, @dls3dml_pnplu, @dlsu3dml_pnplu};   
    method_labels = {'P3P', 'DLS', 'DLSU', 'DLSL', 'DLSLU'};
    
%     method_list = {@p3p_pnplu, @dlsl_pnplu, @dlslu_pnplu, @epnpl_pnplu, @epnpu_pnplu};   
% %     method_list = {@p3p_pnplu, @dls_pnplu, @dlsu_pnplu, @dls3dml_pnplu, @dlsu3dml_pnplu};   
%     method_labels = {'P3P', 'DLSL', 'DLSU', 'EPnPL', 'EPnPU'};
    
    linestyles = {'-', '-', '--', '-', '--', '-'};
    colors = {'r', 'g', 'g', 'b', 'b', 'k'};    
    
    is_ball = false;
    
    noises_3d = [0.01,   0.5,  0.9];
    for noise_it = 1:length(noises_3d)        
        k_2d = 1.0;
        k_3d = 1.0;
        s_3d_min = 0.01;
        s_3d_max = 0.01+noises_3d(noise_it);
        s_2d = 2;
        s_2dl = 2;

        n = 20;
        npt = n;

        T = 1000;

        stats = zeros(2, T, length(method_list));
        p3p_fail_cnt = 0;

        for test_id = 1:T
            %generate rotation and translation
            R = rodrigues(10*randn(3,1));
            t = rand(3, 1);

            %generate 3D points in cam coordinates
            Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[2 10])];
            X = R'*(Xc-repmat(t,1,n));        
            
            Xcs = [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[2 10])];
            Xs = R'*(Xcs-repmat(t,1,n));        
            Xce = [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[2 10])];
            Xe = R'*(Xce-repmat(t,1,n));        

            %generate covariances and perturbed points
            Sigmas3DLines = zeros(6, 6, n);
            Xu = zeros(3, n);
            Xsu = zeros(3, n);
            Xeu = zeros(3, n);
            for i = 1:n                
                
                [R_cov,~,~] = svd(randn(6,6));
                Scale_mat = diag(s_3d_min*ones(6,1) + rand(6,1)*(s_3d_max-s_3d_min));                
                                
                Sigmas3DLines(:,:,i) = R_cov * Scale_mat.^2 * R_cov';                
                noise_vec = R_cov*Scale_mat*randn(6,1);
                
                if is_ball                    
                    scalevec = diag(Scale_mat);
                    alpha = sqrt(mean(scalevec.^2));
                    SavedScale = alpha*eye(6);
                    Sigmas3DLines(:,:,i) = SavedScale*SavedScale;
                    noise_vec = SavedScale*randn(6,1);
                end
                Xsu(:,i) = Xs(:,i) + noise_vec(1:3);
                Xeu(:,i) = Xe(:,i) + noise_vec(4:6);                                
                                R_cov = rodrigues(10*randn(3,1));
                Scale_mat = diag(s_3d_min*ones(3,1) + (diag([1, k_3d, k_3d^2])/k_3d^2*rand(3,1)*(s_3d_max-s_3d_min)));                
                
                Sigmas3D(:,:,i) = R_cov*Scale_mat*Scale_mat*R_cov';
                if is_ball
                    scalevec = diag(Scale_mat);
                    alpha = sqrt(mean(scalevec.^2));
                    SavedScale = alpha*eye(3);
                    Sigmas3D(:,:,i) = SavedScale*SavedScale;
                end
                Xu(:, i) = X(:, i) + R_cov*Scale_mat*randn(3,1);        
            end

            %project perturbed points
            x = R*Xu + repmat(t, 1, n);    
            x = x./repmat(x(3, :), 3, 1);    

            %camera calibration (100 px, fov=90 deg.)
            f = 1000;
            K = [f 0 f/2;
                0 f f/2; 
                0 0 1];            
            iK = inv(K);            
            x = K*x;

            x = x(1:2,:);

            %add projection noise
            xu = zeros(2, n);
            Sigmas2D = zeros(2, 2, n);
            for i = 1:n
                theta = rand(1)*2*pi;
                R_2d = [cos(theta) sin(theta);
                        -sin(theta) cos(theta)];
                Scale_2d = diag(s_2d*diag([1, k_2d])*rand(2,1));
                Sigmas2D(:,:,i) = R_2d*Scale_2d*Scale_2d*R_2d';
                xu(:, i) = x(:, i) + R_2d*Scale_2d*randn(2,1);
            end
            
            xs = R*Xsu + repmat(t, 1, n);    
            xs = xs./repmat(xs(3, :), 3, 1);    
            xs = K*xs;
            xs = xs(1:2,:);
            
            xe = R*Xeu + repmat(t, 1, n);    
            xe = xe./repmat(xe(3, :), 3, 1);    
            xe = K*xe;
            xe = xe(1:2,:);
            
            lu = zeros(3,n);                        
            Sigmas2DLines = zeros(3, 3, n);
            xsus = zeros(2, n);
            xeus = zeros(2, n);
            for i = 1:n                
                theta = rand(1)*2*pi;
                R_2ds = [cos(theta) sin(theta);
                        -sin(theta) cos(theta)];                    
                Scale_2ds = diag(s_2dl*diag([1, k_2d])*rand(2,1));
%                 Scale_2ds = eye(2)*rand;
                theta = rand(1)*2*pi;
                R_2de = [cos(theta) sin(theta);
                        -sin(theta) cos(theta)];                    
                Scale_2de = diag(s_2dl*diag([1, k_2d])*rand(2,1));
%                 Scale_2de = eye(2)*rand;
                xsu = xs(:,i) + R_2ds*Scale_2ds*randn(2,1);
                xeu = xe(:,i) + R_2de*Scale_2de*randn(2,1);
                xsus(:, i) = xsu;
                xeus(:, i) = xeu;
                xsu_h = [xsu;1];
                xeu_h = [xeu;1];
                
                %compute cov directly for a normalized line
                l_h = cross(xsu_h, xeu_h);
                l_h_norm = norm(l_h(1:2));
                l = l_h / l_h_norm;
                
                Sigma_xsu = trace(R_2ds * Scale_2ds.^2 * R_2ds')*eye(2)*0.5;
                Sigma_xeu = trace(R_2de * Scale_2de.^2 * R_2de')*eye(2)*0.5;                
                                
                Sigmas2DLines(:,:,i) = get_line2D_covariance(Sigma_xsu, Sigma_xeu, xsu_h, xeu_h);
                Sigmas2DLinesn(:,:,i) = get_line2D_covariance(Sigma_xsu/(f^2), Sigma_xeu/(f^2), K\xsu_h, K\xeu_h);
                
                ln = K'*l;
                ln = ln/norm(ln(1:2));
                lu(:,i) = ln;
            end
                                    
            %call p3p to get R_est, t_est
            p3p_inds = randsample(n,3);
            X_3 = X(:, p3p_inds);
            xu_3 = xu(:, p3p_inds);
            xuh_3 = iK * [xu_3; ones(1,3)];        
            [R_est t_est] = p3p_interface(X_3, xuh_3);   
            n_sol = size(t_est, 2);
            if n_sol == 0
                p3p_fail_cnt = p3p_fail_cnt  + 1;
                continue
            end
            pose_errors = zeros(n_sol,1);
            for sol_ind = 1:n_sol
                Rc = R_est(:,:,sol_ind);
                tc = t_est(:, sol_ind);
                pose_errors(sol_ind) = norm([Rc',  -Rc'*tc; 0 0 0 1]* [R t; 0 0 0 1] - eye(4));
            end        
            [best_err, best_sol] = min(pose_errors);
            R_est = R_est(:, :, best_sol);
            t_est = t_est(:, best_sol);
            %call PnPU methods        
            for meth_id = 1:length(method_list)                               
                    [Rs, ts] = method_list{meth_id}(X, xu, K, Sigmas3D, Sigmas2D, ...
                        R_est, t_est, Xs, Xe, lu, Sigmas3DLines, Sigmas2DLinesn, xsus, xeus);
%                     [Rs, ts] = method_list{meth_id}(X, xu, K, Sigmas3D, Sigmas2D, ...
%                         R, t, Xs, Xe, lu, Sigmas3DLines, Sigmas2DLinesn, xsus, xeus);

                if size(ts,2) < 1
                    fprintf(['   The solver - ',num2str(meth_id),' - returns no solution! \n']);                                
                    continue;
                end

                %find the best solution            
                error = inf;
                for jjj = 1:size(Rs,3)
                    tempy = cal_pose_err([Rs(:,:,jjj) ts(:,jjj)],[R t]);
                    if sum(tempy) < error
                        y = tempy;
                        error = sum(tempy);
                    end
                end
                stats(:, test_id, meth_id) = y;            
            end
        end
        fprintf('p3p failed in %d cases\n', p3p_fail_cnt);
        disp('Method: median rot/mean rot median trans/mean trans');
        for meth_id = 1:length(method_list)
            r_errs = stats(1, :, meth_id);
            t_errs = stats(2, :, meth_id);
            cnt_ok = length(find(r_errs>0));
            share_ok = int32(cnt_ok / (T-p3p_fail_cnt)*100);
            med_r = median(r_errs);
            mean_r = mean(r_errs);
            std_r = std(r_errs);
            med_t = median(t_errs); 
            mean_t = mean(t_errs);
            std_t = std(t_errs);
            fprintf('%s: %d, %f/%f %f/%f\n', method_labels{meth_id}, share_ok, med_r, mean_r, med_t, mean_t);
            
            report(noise_it, meth_id, :) = [med_r, mean_r, std_r, med_t, mean_t, std_t];
        end        
    end
   save('noise3d_lines_report.mat', 'report', 'method_labels', 'noises_3d');
end