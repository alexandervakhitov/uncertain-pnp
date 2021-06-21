function [err, fail_share, errs2, y_med] = me_kitti_error(fpath, gt_path)
    T = read_kitti_trajectory(fpath, false);    
    T_gt = read_kitti_trajectory(gt_path, false);    
    n = size(T, 3);    
    T_gt = T_gt(:, :, 1:n+1);    
    T_old = [T_gt(:,:,1); 0 0 0 1];    
    for i = 2:n+1
        T_curr = T_gt(:, :, i);
        T_curr_ext = [T_curr; 0 0 0 1];
        T_inc = inv(T_curr_ext)*T_old;
        T_gt(:,:,i-1) = T_inc(1:3,:);
        T_old = T_curr_ext;        
    end
    
    y = zeros(n,2);
    for i = 1:n
        y(i, :) = cal_pose_err(T(:,:,i), T_gt(:,:,i));
    end
    
    y_med = median(y);
    
    ts = reshape(T(:, 4, :), 3, n);
    ts_norms = sqrt(sum(ts.^2, 1));
    inds = find(ts_norms>0);
    ts = ts(:, inds);
    n_inds = length(inds);
    ts_gt = reshape(T_gt(:, 4, inds), 3, n_inds);
    errs2 = sum((ts-ts_gt).^2, 1);
    err = mean(sqrt(errs2));
    fail_share = n_inds/n;
end