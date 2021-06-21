function [R,t] = MLPNP_without_COV(XX, xx, cov)
T = MLPnP(XX, normc(xx));
R = T(1:3,1:3);
t = T(1:3,4);
end

