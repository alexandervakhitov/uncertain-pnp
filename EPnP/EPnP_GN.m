function [R,t,XX0] = EPnP_GN(XX,xx,is_pca)
if nargin < 3
    is_pca = true;
end
Xc = mean(XX, 2);
XX0 = XX;
XX = XX-repmat(Xc, 1, size(XX, 2));
[R,t]= efficient_pnp_gauss(XX.',xx.',diag([1 1 1]), is_pca);
t = t - R*Xc;
return