function Cw = define_control_points_pca(Xw)

Cp1 = mean(Xw, 2);
n = size(Xw, 2);
C = zeros(3,3);
for i = 1:n
    C = C + (Xw(:,i)-Cp1) * (Xw(:, i)-Cp1)';
end
[U,S,V] = svd(C);
Cw = [(V*sqrt(S)/sqrt(n)+repmat(Cp1, 1, 3))';
Cp1'];
