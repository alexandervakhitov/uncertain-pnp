function q = rot2quat(R)
    t = trace(R);
    r = sqrt(1+t);
    w = 0.5*r;
    x = cs(0.5*sqrt(1+R(1,1)-R(2,2)-R(3,3)), R(3,2)-R(2,3));
    y = cs(0.5*sqrt(1-R(1,1)+R(2,2)-R(3,3)), R(1,3)-R(3,1));
    z = cs(0.5*sqrt(1-R(1,1)-R(2,2)+R(3,3)), R(2,1)-R(1,2));
    q = [w;x;y;z];
end
function cs = cs(x,y)
    cs = sign(y)*abs(x);
end