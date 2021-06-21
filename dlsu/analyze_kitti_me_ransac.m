function analyze_kitti_me_ransac
    addpath('/home/alexander/materials/line_descriptor/code/lld/matlab');

    epnp_path = '/home/alexander/materials/pnp3d/ransac/epnp_flat/0/_epnp_traj.txt';
    epnp_inl_path = '/home/alexander/materials/pnp3d/ransac/epnp_flat/0/_epnp_inliers.txt';
    
    epnpu_label = 'epnpu';
    
    epnpu_path = ['/home/alexander/materials/pnp3d/ransac/epnp_flat/0/_' epnpu_label '_traj.txt'];
    epnpu_inl_path = ['/home/alexander/materials/pnp3d/ransac/epnp_flat/0/_' epnpu_label '_inliers.txt'];
    gt_path = '/home/alexander/materials/sego/kitti_odometry/dataset/poses/00.txt';
    n = 2;
    errors = zeros(n,1);
    fail = zeros(n,1);            
    
    [err, good_share, err2_epnpu, y_med_epnpu] = me_kitti_error(epnpu_path, gt_path);
    inliers_epnpu = dlmread(epnpu_inl_path);
    errors(2) = err;
    fail(2) = 1-good_share;
    
    [err, good_share, err2_epnp, y_med_epnp] = me_kitti_error(epnp_path, gt_path);
    inliers_epnp = dlmread(epnp_inl_path);
    errors(1) = err;
    fail(1) = 1-good_share;
    
    Y = [y_med_epnp;y_med_epnpu];
    I = [mean(inliers_epnp, 1); mean(inliers_epnpu, 1)];
    [errors fail Y I]
    N = length(err2_epnp);
    plot(1:N, err2_epnp, 'b', 1:N, err2_epnpu, 'r');
end