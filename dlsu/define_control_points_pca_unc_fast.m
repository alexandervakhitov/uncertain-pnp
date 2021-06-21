function Cw = define_control_points_pca_unc_fast(Xw, sigmas3D)
n = size(Xw, 2);
Cp1 = zeros(3,1);
sigmas3D = sigmas3D + 1e-10;
S_sum = 0;
for i = 1:n
    Sci = 1.0/(sigmas3D(i));
    Cp1 = Cp1 + Sci*Xw(:, i);
    S_sum = S_sum + Sci;
end
Si = 1.0/S_sum;
Cp1 = Si*Cp1;

C = zeros(3,3);
for i = 1:n       
    Xmc = (Cp1-Xw(:,i));
    C = C + Xmc * Xmc' / (sigmas3D(i) - Si);    
end
[U, S, V] = svd(C);
Cw = [(V*sqrt(S)/sqrt(n)+repmat(Cp1, 1, 3))';
Cp1'];
