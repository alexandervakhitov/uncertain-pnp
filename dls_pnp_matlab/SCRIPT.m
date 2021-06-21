% This script file demonstrates how to call the dls pnp solver

clear; clc; close all;

% number of points
n = 10;

% camera pixel noise sigma in normalized pixel coordinates (i.e., in this case we have sigma = 1 pixel with a focal length of 800 pixels)
sigma = 1 * 1/800;

% real world simulator (rws)
% generate n random points (p), and their projections on the image (z), the
% global-to-camera transformation is the rotation matrix C (global expressed w.r.t. camera), and the
% translation t (origin of global w.r.t. camera)
[C, t, p, z] = rws(n, sigma);

% call the dls_pnp solver
% returns the estimated rotation matrix (Cm), and the estimated translation
% vector (tm) along with the cost of each solution, and a flag incase there
% was an error

%[Cm, tm, cost, flag] = dls_pnp(p, z);
[Cm, tm, cost, flag] = robust_dls_pnp(p, z); % version that calls dls 3 times, to avoid Cayley singularity

for i = 1:size(tm,2)
    % compute the error w.r.t. the true, and display the solution
    [da, dt] = compute_error(C, t, Cm(:,:,i), tm(:,i));
    % Error is
    fprintf('DLS PnP sol %d: orientation error = %.3f deg, translation error = %.3f m, cost = %f\n', i, da * 180/pi, dt, cost(i))
    
%    plot_reprojection(z, p, Cm(:,:,i), tm(:,i));
end

