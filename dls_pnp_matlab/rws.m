function [C, t, p, z] = rws(N, sigma)
% this function generates a random camera pose, along with N random points,
% and also the perspective projections of those points.

% Generate a random global-to-camera rotation. This is the orientation of 
% the global frame expressed in the camera frame of refence.
angle = 15;

C = rotx(angle * randn * pi/180 ) * roty( angle * randn * pi/180 );

% Generate a random global-to-camera translation. This is the origin of the
% global frame expressed in the camera frame.
t = randn(3,1);

% Create random 3D points within a 45 deg FOV (vertical and horizontal) of
% the camera. The points are between 0.5 and 5.5 meters from the camera.
Psens = zeros(3,N);
theta = (rand(N,1)*45 - 22.5) * pi/180;
phi = (rand(N,1)*45 - 22.5) * pi/180;
for i = 1:N
    psens_unit = rotx(theta(i)) * roty(phi(i)) * [0;0;1];
    alpha = rand * 5 + 0.5;
    Psens(:,i) = alpha * psens_unit;
end

% Express the points in the global frame of reference
p = C' *(Psens - repmat(t,1,N));

% Construct the vector of perspective projections (i.e., image
% measurements) of the points,
z = zeros(2,N);
for i = 1:N
    
    % create an instance of 2x1 pixel noise
    noise = sigma * randn(2,1);
    
    % You can uncomment the following lines in order to limit the noise to +/-
    % 3 sigma
    %
    %
    %     if abs(noise(1)) > 3 * sigma
    %         noise(1) = sign(noise(1)) * 3 * sigma;
    %     end
    %     if abs(noise(2)) > 3 * sigma
    %         noise(2) = sign(noise(2)) * 3 * sigma;
    %     end
    
    % Create the image measurement using the standard pinhole camera model
    z(:,i) = [ Psens(1,i) / Psens(3,i) ; Psens(2,i) / Psens(3,i)] + noise;
end

end

function r = rotx(t)
%rotx: rotation around the x-axis

ct = cos(t);
st = sin(t);
r =    [1	0	0;
    0	ct	-st;
    0	st	ct];
end

function r = roty(t)
% roty: rotation about y-axis
ct = cos(t);
st = sin(t);
r =    [ct	0	st;
    0	1	0;
    -st	0	ct];
end
