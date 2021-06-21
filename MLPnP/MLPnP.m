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

% 28.06.2016 by Steffen Urban
% if you use this file it would be neat to cite our paper:
% @INPROCEEDINGS {mlpnp2016,
%    author    = "Urban, S.; Leitloff, J.; Hinz, S.",
%    title     = "MLPNP - A REAL-TIME MAXIMUM LIKELIHOOD SOLUTION TO THE PERSPECTIVE-N-POINT PROBLEM.",
%    booktitle = "ISPRS Annals of Photogrammetry, Remote Sensing \& Spatial Information Sciences",
%    year      = "2016",
%    volume    = "3",
%    pages     = "131-138"}
%% MLPnP - Maximum Likelihood Perspective-N-Point

% input: 1. points3D - a 3xN matrix of N 3D points in the object coordinate system
%        2. v        - a 3xN matrix of N bearing vectors (camera rays)
%                      ||v|| = 1
%        3. cov      - if covariance information of bearing vectors if
%                      available then cov is a 9xN matrix. 
%                      e.g. it can be computed from image plane variances 
%                      sigma_x and sigma_y (in case of perspective cameras): 
%                      cov = K\diag([sigma_x sigma_y 0])/K'
%                      cov = reshape(cov,9,1)
%                      here K\ and /K' are the Jacobians of the image to
%                      bearing vector transformation (inverse calibration
%                      matrix K. Details in the paper.
% output: 1. T          - 4x4 transformation matrix [R T;0 0 0 1]
%         2. statistics - contains statistics after GN refinement, see
%                         optim_GN.m for details

function [T, r, s, Kll] = MLPnP(points3D, v, cov, is_ref)
if nargin < 4
   is_ref = true; 
end

use_cov = 1;
% if cov is not given don't use it
if nargin < 3
    use_cov = 0;
end

nrPts = size(points3D,2);
% matrix of null space vectors r and s
r = zeros(3,nrPts);
s = zeros(3,nrPts);
cov_reduced = zeros(2,2,nrPts);

% test planarity, only works well if the scene is really planar
% quasi-planar won't work very well
S = points3D*points3D';
[eigRot,~] = eig(S);
planar = 0;
% create full design matrix
A = zeros(nrPts,12);
if (rank(S) == 2)
    planar = 1;
    points3D1 = eigRot'*(points3D);
    points3Dn = [points3D1;ones(1,nrPts)];
    % create reduced design matrix
    A = zeros(nrPts,9);
else
    points3Dn = [points3D;ones(1,nrPts)];
end

% compute null spaces of bearing vector v: null(v')
for i=1:nrPts
    null_2d = null(v(1:3,i)');
    r(:,i) = null_2d(:,1);
    s(:,i) = null_2d(:,2);
    if use_cov
        tmp = reshape(cov(:,i),3,3);
        cov_reduced(:,:,i) = (null_2d'*tmp*null_2d)^-1;
    end
end

% stochastic model
Kll = eye(2*nrPts,2*nrPts);
% if (normalize)
%     points3Dn = normc(points3Dn);
% end
if planar % build reduces system
    for i=1:nrPts
        if (use_cov)
            Kll(2*i-1:2*i,2*i-1:2*i) = cov_reduced(:,:,i);
        end
        % r12
        A (2*i-1,1) = r(1,i)*points3Dn(2,i);
        A (2*i,1) = s(1,i)*points3Dn(2,i);
        % r13
        A (2*i-1,2) = r(1,i)*points3Dn(3,i);
        A (2*i,2) = s(1,i)*points3Dn(3,i); 
        % r22
        A (2*i-1,3) = r(2,i)*points3Dn(2,i);
        A (2*i,3) = s(2,i)*points3Dn(2,i);
        % r23
        A (2*i-1,4) = r(2,i)*points3Dn(3,i);
        A (2*i,4) = s(2,i)*points3Dn(3,i);
        % r31
        A (2*i-1,5) = r(3,i)*points3Dn(2,i);
        A (2*i,5) = s(3,i)*points3Dn(2,i);
        % r32
        A (2*i-1,6) = r(3,i)*points3Dn(3,i);
        A (2*i,6) = s(3,i)*points3Dn(3,i);   
        % t1
        A (2*i-1,7) = r(1,i);
        A (2*i,7)   = s(1,i);
        % t2
        A (2*i-1,8) = r(2,i);
        A (2*i,8)   = s(2,i);
        % t3
        A (2*i-1,9) = r(3,i);
        A (2*i,9)   = s(3,i);

    end 
else % build full system
    for i=1:nrPts
       if (use_cov)
            Kll(2*i-1:2*i,2*i-1:2*i) = cov_reduced(:,:,i);
       end
       % r11
       A (2*i-1,1) = r(1,i)*points3Dn(1,i);
       A (2*i,1) = s(1,i)*points3Dn(1,i);
       % r12
       A (2*i-1,2) = r(1,i)*points3Dn(2,i);
       A (2*i,2) = s(1,i)*points3Dn(2,i); 
       % r13
       A (2*i-1,3) = r(1,i)*points3Dn(3,i);
       A (2*i,3) = s(1,i)*points3Dn(3,i);
       % r21
       A (2*i-1,4) = r(2,i)*points3Dn(1,i);
       A (2*i,4) = s(2,i)*points3Dn(1,i);
       % r22
       A (2*i-1,5) = r(2,i)*points3Dn(2,i);
       A (2*i,5) = s(2,i)*points3Dn(2,i);
       % r23
       A (2*i-1,6) = r(2,i)*points3Dn(3,i);
       A (2*i,6) = s(2,i)*points3Dn(3,i);
       % r31
       A (2*i-1,7) = r(3,i)*points3Dn(1,i);
       A (2*i,7) = s(3,i)*points3Dn(1,i);
       % r32
       A (2*i-1,8) = r(3,i)*points3Dn(2,i);
       A (2*i,8) = s(3,i)*points3Dn(2,i);
       % r33
       A (2*i-1,9) = r(3,i)*points3Dn(3,i);
       A (2*i,9) = s(3,i)*points3Dn(3,i);    
       % t1
       A (2*i-1,10) = r(1,i);
       A (2*i,10)   = s(1,i);
       % t2
       A (2*i-1,11) = r(2,i);
       A (2*i,11)   = s(2,i);
       % t3
       A (2*i-1,12) = r(3,i);
       A (2*i,12)   = s(3,i); 
    end
end

% do least squares AtPAx=0
b = A'*A;    
[~,~,v1] = svd(b);

if planar
    tout1  = v1(7:9,end);
    P = zeros(3,3);
    P(:,2:3) = reshape(v1(1:6,end),2,3)';
    scalefact = sqrt(abs(norm(P(:,2))*norm(P(:,3))));
    P(:,1) = cross(P(:,2),P(:,3));
    P = P';
     %SVD to find the best rotation matrix in the Frobenius sense
    [U2,~,V2] = svd(P(1:3,1:3));
    R = U2*V2'; 
    if det(R) < 0
        R = -1*R;
    end
    % rotate solution back (see paper)
    R = eigRot*R;
    % recover translation
    tout = (tout1./scalefact);
    R = -R';
    
    R1 = [R(:,1) R(:,2) R(:,3)];
    R2 = [-R(:,1) -R(:,2) R(:,3)];
    Ts = zeros(4,4,4);
    Ts(:,:,1) = [R1 tout;0 0 0 1];
    Ts(:,:,2) = [R1 -tout;0 0 0 1];
    Ts(:,:,3) = [R2 tout;0 0 0 1];
    Ts(:,:,4) = [R2 -tout;0 0 0 1];
    % find the best solution with 6 correspondences
    diff1 = zeros(4,1);
    for te = 1:6
        for ba = 1:4
            testres1 = Ts(:,:,ba)*[points3D(:,te);1];
            testres11 = normc(testres1(1:3));
            diff1(ba) = diff1(ba) + (1-dot(testres11,v(:,te)));
        end
    end
    [~,idx] = min(diff1);
    T = Ts(:,:,idx);
else
    tout1  = v1(10:12,end);   
    P = reshape(v1(1:9,end),3,3);
    scalefact = (abs(norm(P(:,1))*norm(P(:,2))*norm(P(:,3))))^(1/3);
    %SVD to find the best rotation matrix in the Frobenius sense
    [U2,~,V2] = svd(P(1:3,1:3)); 
    R = U2*V2';
    if det(R) < 0
        R = -1*R;
    end
    % recover translation
    tout = R*(tout1./scalefact);
    T1 = [R tout;0 0 0 1]^-1;
    T2 = [R -tout;0 0 0 1]^-1;
    diff1 = 0;
    diff2 = 0;
    % find the best solution with 6 correspondences
    for te = 1:6
        testres1 = T1*[points3D(:,te);1];
        testres2 = T2*[points3D(:,te);1];
        testres1 = normc(testres1(1:3));
        testres2 = normc(testres2(1:3));

        diff1 = diff1+(1-dot(testres1,v(:,te)));
        diff2 = diff2+(1-dot(testres2,v(:,te)));
    end
    if diff1 < diff2
        T = T1(1:3,1:4);
    else
        T = T2(1:3,1:4);  
    end
end

if is_ref
    optimFlags.epsP  = 1e-6;
    optimFlags.epsF  = 1e-6;
    optimFlags.maxit = 5;
    optimFlags.tau   = 1e-4;
    [T, statistics] = optim_MLPnP_GN(T, points3D, r, s, Kll, optimFlags);
end
end

