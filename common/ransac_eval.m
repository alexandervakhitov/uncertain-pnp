function  [inliers, R, T] = ransac_eval(K, tR, tT, X, x, t, scales_2d)

R = [];
T = [];
inliers = [];
if size(scales_2d, 1) == 0
    scales_2d = ones(1, size(X, 2));
end

try
%K = [f 0 0;0 f 0; 0 0 1];
for i = 1:size(tR,3)
    newX  = K * (tR(:,:,i) * X + tT(:,i) * ones(1,size(X,2)));
    error = sqrt(sum((x - newX(1:2,:)./repmat(newX(3,:),2,1)).^2)) ./ scales_2d;
    tinliers = find((error < t) & (newX(3, :) > 0));
    if length(tinliers) > length(inliers)
        inliers = tinliers;
        R = tR(:,:,i);
        T = tT(:,i);
    end
end
catch
    R = [];
    T = [];
    inliers = [];
end

end
