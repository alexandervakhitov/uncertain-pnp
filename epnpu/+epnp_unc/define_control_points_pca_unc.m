function Cw = define_control_points_pca_unc(Xw, Sigmas3D)
n = size(Xw, 2);
%Cp1 = mean(Xw, 2);
Cp1 = zeros(3,1);
S_sum = zeros(3,3);
for i = 1:n
    Sci = inv(Sigmas3D(:,:,i)+1e-9*eye(3));
    Cp1 = Cp1 + Sci*Xw(:, i);
    S_sum = S_sum + Sci;
end
Si = inv(S_sum);
Cp1 = Si*Cp1;

C = zeros(3,3);
for i = 1:n
    Sci = sqrt(inv(Sigmas3D(:,:,i)+1e-9*eye(3) - Si));
    Xmc = Sci*(Cp1-Xw(:,i));
    C = C + Xmc * Xmc';
end
[U,S,V] = svd(C);
Cw = [(V*sqrt(S)/sqrt(n)+repmat(Cp1, 1, 3))';
Cp1'];
