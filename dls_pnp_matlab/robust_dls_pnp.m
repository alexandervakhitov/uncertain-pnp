function [C_est, t_est, cost, flag] = robust_dls_pnp(p, z)

C_est = [];
t_est = [];
cost = 1e10;
flag = 0;

C_temp1 = zeros(3,3,0); t_temp1 = zeros(3,0);
C_temp2 = zeros(3,3,0); t_temp2 = zeros(3,0);

R = cat(3, rotx(pi/2), roty(pi/2), rotz(pi/2));
t = mean(p,2);

cost = inf;

for i = 1:3
    % Make a random rotation
    pp = R(:,:,i) * (p - repmat(t, 1, size(p,2)));
    
    timer_solver = tic;
    [C_est_i, t_est_i, cost_i, flag_i] = dls_pnp(pp, z);
    time_solve = toc(timer_solver);
%     fprintf('dls solver time %f\n', time_solve);
    for j = 1:length(cost_i)
        t_est_i(:,j) = t_est_i(:,j) - C_est_i(:,:,j) * R(:,:,i) * t;
        C_est_i(:,:,j) = C_est_i(:,:,j) * R(:,:,i);
    end
    
%     if min(cost_i) < min(cost)
%         C_est = C_est_i;
%         t_est = t_est_i;
%         cost = cost_i;
%         flag = flag_i;
%     end
    min_cost_i = min(cost_i);
    index = find(cost_i < min_cost_i*2);
    C_temp1 = cat(3,C_temp1, C_est_i(:,:,index));
    t_temp1 = [t_temp1 t_est_i(:,index)];    
    
    C_temp2 = cat(3,C_temp2, C_est_i);
    t_temp2 = [t_temp2 t_est_i];
end

if size(p,2) <= 6
    C_est = C_temp1;
    t_est = t_temp1;
    
else
    %since there are three different objectives in DLS+++,
    %we should use a single objective to select the final solution.
    %here, we use the reprojection error.

    for i = 1:size(C_temp2,3)
        proj = C_temp2(:,:,i)*p + t_temp2(:,i)*ones(1,size(p,2));
        proj = proj(1:2,:)./repmat(proj(3,:),2,1);
        cost_i = norm(z-proj,'fro');
        if cost_i < cost
            C_est = C_temp2(:,:,i);
            t_est = t_temp2(:,i);
            cost = cost_i;
        end
    end
end
end

function r = rotx(t)
ct = cos(t);
st = sin(t);
r =    [1	0	0;
    0	ct	-st;
    0	st	ct];
end

function r = roty(t)
% roty: rotation about y-axi-
ct = cos(t);
st = sin(t);
r =    [ct	0	st;
    0	1	0;
    -st	0	ct];
end

function r = rotz(t)
% rotz: rotation about z-axis

ct = cos(t);
st = sin(t);
r = [ct	-st	0
    st	ct	0
    0	0	1];

end
