%     Steffen Urban email: urbste@googlemail.com
%     Copyright (C) 2016  Steffen Urban
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% 28.06.2016 Steffen Urban
function [err,J] = residualsAndJacobianLsq(x, r, s, points3D, ICov)

nrPts = size(points3D,2);

err = zeros(2*nrPts,1);
J = zeros(2*nrPts,6);
R = Rodrigues2(x(1:3));
t = x(4:6);

res1 = R*points3D+repmat(t,1,nrPts);
normres = normc(res1(1:3,:));

for i=1:size(r,2)
    err(2*i-1,1) = r(:,i)'*normres(:,i);
    err(2*i,1) = s(:,i)'*normres(:,i);
    J(2*i-1:2*i,1:6) = jacobians_Rodrigues(points3D(1,i),points3D(2,i),points3D(3,i),...
        r(1,i),r(2,i),r(3,i),s(1,i),s(2,i),s(3,i),x(4),x(5),x(6),x(1),x(2),x(3));
end
    J = ICov * J;
    err = ICov * err;
end

