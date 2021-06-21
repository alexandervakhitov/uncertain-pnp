function data = read_kitti_trajectory(fpath, is_inv)
    data = dlmread(fpath);
    n = size(data, 1);
    data = reshape(data, n, 4, 3);
    data = permute(data, [3, 2, 1]);  
    if (is_inv)
        n = size(data, 3);
        for i = 1:n
            T0 = data(:, :, i);
            T1 = [T0(:, 1:3)', -T0(:, 1:3)'*T0(:, 4)];
            data(:, :, i) = T1;
        end
    end
end