function pnpu_test_noise2d
    addpath('../');
    addpath('../dls_pnp_matlab');
    addpath('../p3p_code_final');
    addpath('../EPnP_matlab/EPnP');
    addpath('/home/alexander/materials/pnp3d/MLPnP_matlab_toolbox/MLPnP');
    addpath('/home/alexander/materials/pnp3d/upnp/opengv/matlab');

    method_list = {@p3p_pnpu, @dls_pnpu, @dlsu_pnpu, @epnp_pnpu, @epnpu_pnpu, @mlpnp_pnpu};   %, @upnp_pnpu
    method_labels = {'P3P', 'DLS', 'DLSU', 'EPnP', 'EPnPU', 'MLPnP', 'UPnP'};
    linestyles = {'-', '-', '--', '-', '--'};
    colors = {'r', 'g', 'g', 'b', 'b'};    
    
    noises_2d = 2*[0.01, 1.0, 3.0, 5.0, 7.0, 9.0];
    for noise_it = 1:length(noises_2d)

        k_3d = sqrt(10.0);
        k_2d = 1.0;%sqrt(10.0);
        s_3d_min = 0.01;
        s_3d_max = 0.01;
        s_2d = noises_2d(noise_it);

        n = 20;

        T = 1000;

        stats = zeros(2, T, length(method_list));
        p3p_fail_cnt = 0;

        for test_id = 1:T
            %generate rotation and translation
            R = rodrigues(10*randn(3,1));
            t = rand(3, 1);

            %generate 3D points in cam coordinates
            Xc = randn(3, n);
            Xc(3,:) = rand(1, n);    
            Xc(3,:) = 9*Xc(3,:)+1;
            X = R'*(Xc-repmat(t,1,n));        

            %generate covariances and perturbed points
            Sigmas3D = zeros(3, 3, n);
            Xu = zeros(3, n);
            for i = 1:n
                R_cov = rodrigues(10*randn(3,1));
                Scale_mat = diag(s_3d_min*ones(3,1) + (diag([1, k_3d, k_3d^2])/k_3d^2*rand(3,1)*(s_3d_max-s_3d_min)));
                Sigmas3D(:,:,i) = R_cov*Scale_mat*Scale_mat*R_cov';
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
                [Rs, ts] = method_list{meth_id}(X, xu, K, Sigmas3D, Sigmas2D, R_est, t_est);
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
            med_t = median(t_errs); 
            mean_t = mean(t_errs);
            std_r = std(r_errs);
            std_t = std(t_errs);
            fprintf('%s: %d, %f/%f %f/%f\n', method_labels{meth_id}, share_ok, med_r, mean_r, med_t, mean_t);
            
            report(noise_it, meth_id, :) = [med_r, mean_r, std_r, med_t, mean_t, std_t];
        end        
    end
   save('noise2d_report.mat', 'report', 'method_labels', 'noises_2d');
end