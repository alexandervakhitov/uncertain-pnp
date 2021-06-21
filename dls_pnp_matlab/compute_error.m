function [da, dt] = compute_error(C, t, Cm, tm)
% compute the error quaternion btw. the true and the estimated solutions
% (using JPL definition of quaternions)
q_del = rot2quat(C' * Cm);

% compute the tilt angle error
da = norm(q_del(1:3) * 2);

% compute the position error
dt = norm(t - tm);
end

function q = rot2quat(R)
% converts a rotational matrix to a unit quaternion, according to JPL
% procedure (Breckenridge Memo)

T = trace(R);

[dummy maxpivot] = max([R(1,1) R(2,2) R(3,3) T]); %#ok<ASGLU>

switch maxpivot
    case 1
        q(1) = sqrt((1+2*R(1,1)-T)/4);
        q(2:4) = 1/(4*q(1)) * [R(1,2)+R(2,1);
            R(1,3)+R(3,1);
            R(2,3)-R(3,2) ];
        
    case 2
        q(2) = sqrt((1+2*R(2,2)-T)/4);
        q([1 3 4]) = 1/(4*q(2)) * [R(1,2)+R(2,1);
            R(2,3)+R(3,2);
            R(3,1)-R(1,3) ];
        
    case 3
        q(3) = sqrt((1+2*R(3,3)-T)/4);
        q([1 2 4]) = 1/(4*q(3)) * [R(1,3)+R(3,1);
            R(2,3)+R(3,2);
            R(1,2)-R(2,1) ];
        
    case 4
        q(4) = sqrt((1+T)/4);
        q(1:3) = 1/(4*q(4)) * [R(2,3)-R(3,2);
            R(3,1)-R(1,3);
            R(1,2)-R(2,1) ];
        
end % switch

% make column vector
q = q(:);

% 4th element is always positive
if q(4)<0
    q = -q;
end

% quaternion normalization
q = q/sqrt(q'*q);
end