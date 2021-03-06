function dsq=define_distances_btw_control_points_pca(Cw)

c1 = Cw(1,:);
c2 = Cw(2,:);
c3 = Cw(3,:);
c4 = Cw(4,:);

d12=(c1(1)-c2(1))^2 + (c1(2)-c2(2))^2 + (c1(3)-c2(3))^2;
d13=(c1(1)-c3(1))^2 + (c1(2)-c3(2))^2 + (c1(3)-c3(3))^2;
d14=(c1(1)-c4(1))^2 + (c1(2)-c4(2))^2 + (c1(3)-c4(3))^2;
d23=(c2(1)-c3(1))^2 + (c2(2)-c3(2))^2 + (c2(3)-c3(3))^2;
d24=(c2(1)-c4(1))^2 + (c2(2)-c4(2))^2 + (c2(3)-c4(3))^2;
d34=(c3(1)-c4(1))^2 + (c3(2)-c4(2))^2 + (c3(3)-c4(3))^2;

dsq=[d12,d13,d14,d23,d24,d34]';