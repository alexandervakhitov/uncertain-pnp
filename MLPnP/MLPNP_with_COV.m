function [R,t, r, s, Kll] = MLPNP_with_COV(XX,xx,cov,is_full)
if size(xx, 1) == 2
   homx = [xx; ones(1,size(xx,2))];
   xx = normc(homx); 
   sigmas_2d = compute_traces(cov);
   Evv = zeros(3, 3, size(xx, 2));
   cov = zeros(9, size(xx, 2));
   for id = 1:size(xx, 2)
       cov_proj = diag([sigmas_2d(id) sigmas_2d(id) 0]);          
       J = (eye(3)-(xx(:,id)*xx(:,id)')/(xx(:,id)'*xx(:,id)))/norm(homx(:,id));
       Evv(:,:,id) = J*cov_proj*J';
       cov(:,id) = reshape(Evv(:,:,id),9,1);
   end
end
if nargin < 4
    [T, r, s, Kll] = MLPnP(XX, normc(xx), cov);
else
    [T, r, s, Kll] = MLPnP(XX, normc(xx), cov, is_full);
end
R = T(1:3,1:3);
t = T(1:3,4);
end

