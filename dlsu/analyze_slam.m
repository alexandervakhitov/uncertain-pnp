function analyze_slam
    addpath('/home/alexander/materials/line_descriptor/code/lld/matlab');

    folder_path = '/home/alexander/materials/pnp3d/segoexp/ORB_SLAM2_SEGO/RELOC';
    gt_path = '/home/alexander/materials/sego/kitti_odometry/dataset/poses/00.txt';
    n = 5;
    errors = zeros(n,1);
    fail = zeros(n,1);
    for i = 0:n-1
        fpath = [folder_path '/traj_00_' num2str(i) '.txt'];
        [err, good_share] = slam_kitti_error(fpath, gt_path);
        errors(i+1) = err;
        fails(i+1) = 1-good_share;
    end
    median(errors)
    fprintf('%d-%d-%d\n', min(fails), median(fails), max(fails));
    
end