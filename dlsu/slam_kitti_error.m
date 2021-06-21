function [err, fail_share] = slam_kitti_error(fpath, gt_path)
    T = read_kitti_trajectory(fpath, true);    
    T_gt = read_kitti_trajectory(gt_path, false);
    T_gt = T_gt(:, :, 2:2:end);
    n = size(T, 3);
    ts = reshape(T(:, 4, :), 3, n);
    ts_norms = sqrt(sum(ts.^2, 1));
    inds = find(ts_norms>0);
    ts = ts(:, inds);
    n_inds = length(inds);
    ts_gt = reshape(T_gt(:, 4, inds), 3, n_inds);
    err = calc_rmse(ts, ts_gt, true);    
    fail_share = n_inds/n;
end