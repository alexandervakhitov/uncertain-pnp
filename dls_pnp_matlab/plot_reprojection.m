function [] = plot_reprojection(z, p, C, t)


figure;
plot(z(1,:), z(2,:), 'rx'); hold on

cam_pts = C * p + repmat(t, 1, length(p));

im_pts = cam_pts(1:2,:) ./ repmat(cam_pts(3,:), 2,1);

plot(im_pts(1,:), im_pts(2,:), 'bx')