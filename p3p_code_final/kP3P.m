function [R T] = kP3P(Pts, impts )

impts = [impts; ones(1,size(impts,2))];
impts = impts./repmat(sqrt(sum(impts.^2,1)),3,1);

poses = p3p(Pts,impts);

j = 1;
R = [];
T = [];
for i = 1:size(poses,2)/4
    tR = poses(:,4*(i-1)+2:4*i)';
    tT = -tR*poses(:,4*(i-1)+1);
 
    proj = tR*Pts + repmat(tT,1,3);
    if min(proj(3,:)) > 0
        R(:,:,j) = tR;
        T(:,j)   = tT;
        j = j + 1;
    end
end
end
