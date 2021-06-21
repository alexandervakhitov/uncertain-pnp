% Copyright (c) 2011, Laurent Kneip, ETH Zurich
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of ETH Zurich nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL ETH ZURICH BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
% p3p.m
%
%
%      Author: Laurent Kneip
% Description: Compute the absolute pose of a camera using three 3D-to-2D correspondences
%   Reference: A Novel Parametrization of the P3P-Problem for a Direct Computation of
%              Absolute Camera Position and Orientation
%
%       Input: worldPoints: 3x3 matrix with corresponding 3D world points (each column is a point)
%              imageVectors: 3x3 matrix with UNITARY feature vectors (each column is a vector)
%      Output: poses: 3x16 matrix that will contain the solutions
%                     form: [ 3x1 position(solution1) 3x3 orientation(solution1) 3x1 position(solution2) 3x3 orientation(solution2) ... ]
%                     the obtained orientation matrices are defined as transforming points from the cam to the world frame

function poses = p3p( worldPoints, imageVectors )

    % Initialization of the solution matrix, extraction of world points and feature vectors

    poses = zeros(3,4*4);

    P1 = worldPoints(:,1);
    P2 = worldPoints(:,2);
    P3 = worldPoints(:,3);
    
    % Verification that world points are not colinear

    vector1 = P2 - P1;
    vector2 = P3 - P1;

    if norm(cross(vector1,vector2)) == 0
        return
    end

    % Creation of intermediate camera frame

    f1 = imageVectors(:,1);
    f2 = imageVectors(:,2);
    f3 = imageVectors(:,3);

    e1 = f1;
    e3 = cross(f1,f2);
    norm_e3 = sqrt(sum(e3.*e3));
    e3 = e3 ./ repmat(norm_e3,3,1);
    e2 = cross(e3,e1);

    T = [ e1'; e2'; e3' ];

    f3 = T*f3;
    
    % Reinforce that f3[2] > 0 for having theta in [0;pi]
    
    if ( f3(3,1) > 0 )
        
        f1 = imageVectors(:,2);
        f2 = imageVectors(:,1);
        f3 = imageVectors(:,3);

        e1 = f1;
        e3 = cross(f1,f2);
        norm_e3 = sqrt(sum(e3.*e3));
        e3 = e3 ./ repmat(norm_e3,3,1);
        e2 = cross(e3,e1);

        T = [ e1'; e2'; e3' ];

        f3 = T*f3;
        
        P1 = worldPoints(:,2);
        P2 = worldPoints(:,1);
        P3 = worldPoints(:,3);
        
    end
    
    % Creation of intermediate world frame
    
    n1 = P2-P1;
    norm_n1 = sqrt(sum(n1.*n1));
    n1 = n1 ./ repmat(norm_n1,3,1);
    n3 = cross(n1,(P3-P1));
    norm_n3 = sqrt(sum(n3.*n3));
    n3 = n3 ./ repmat(norm_n3,3,1);
    n2 = cross(n3,n1);

    N = [ n1'; n2'; n3' ];
    
    % Extraction of known parameters
    
    P3 = N*(P3-P1);
    
    d_12 = sqrt(sum((P2 - P1).^2));
    f_1 = f3(1,1)/f3(3,1);
    f_2 = f3(2,1)/f3(3,1);
    p_1 = P3(1,1);
    p_2 = P3(2,1);
    
    cos_beta = f1'*f2;
    b = 1/(1-cos_beta^2) - 1;
    
    if cos_beta < 0
        b = -sqrt(b);
    else
        b = sqrt(b);
    end
    
    % Definition of temporary variables for avoiding multiple computation
    
    f_1_pw2 = f_1^2;
    f_2_pw2 = f_2^2;
    p_1_pw2 = p_1^2;
    p_1_pw3 = p_1_pw2 * p_1;
    p_1_pw4 = p_1_pw3 * p_1;
    p_2_pw2 = p_2^2;
    p_2_pw3 = p_2_pw2 * p_2;
    p_2_pw4 = p_2_pw3 * p_2;
    d_12_pw2 = d_12^2;
    b_pw2 = b^2;
    
    % Computation of factors of 4th degree polynomial
    
    factor_4 = -f_2_pw2*p_2_pw4 ...
               -p_2_pw4*f_1_pw2 ...
               -p_2_pw4;

    factor_3 = 2*p_2_pw3*d_12*b ...
               +2*f_2_pw2*p_2_pw3*d_12*b ...
               -2*f_2*p_2_pw3*f_1*d_12;

    factor_2 = -f_2_pw2*p_2_pw2*p_1_pw2 ...
               -f_2_pw2*p_2_pw2*d_12_pw2*b_pw2 ...
               -f_2_pw2*p_2_pw2*d_12_pw2 ...
               +f_2_pw2*p_2_pw4 ...
               +p_2_pw4*f_1_pw2 ...
               +2*p_1*p_2_pw2*d_12 ...
               +2*f_1*f_2*p_1*p_2_pw2*d_12*b ...
               -p_2_pw2*p_1_pw2*f_1_pw2 ...
               +2*p_1*p_2_pw2*f_2_pw2*d_12 ...
               -p_2_pw2*d_12_pw2*b_pw2 ...
               -2*p_1_pw2*p_2_pw2;

    factor_1 = 2*p_1_pw2*p_2*d_12*b ...
               +2*f_2*p_2_pw3*f_1*d_12 ...
               -2*f_2_pw2*p_2_pw3*d_12*b ...
               -2*p_1*p_2*d_12_pw2*b;

    factor_0 = -2*f_2*p_2_pw2*f_1*p_1*d_12*b ...
               +f_2_pw2*p_2_pw2*d_12_pw2 ...
               +2*p_1_pw3*d_12 ...
               -p_1_pw2*d_12_pw2 ...
               +f_2_pw2*p_2_pw2*p_1_pw2 ...
               -p_1_pw4 ...
               -2*f_2_pw2*p_2_pw2*p_1*d_12 ...
               +p_2_pw2*f_1_pw2*p_1_pw2 ...
               +f_2_pw2*p_2_pw2*d_12_pw2*b_pw2;

    % Computation of roots

    x=solveQuartic([factor_4 factor_3 factor_2 factor_1 factor_0]);
    
    % Backsubstitution of each solution
    
    for i=1:4

        cot_alpha = (-f_1*p_1/f_2-real(x(i))*p_2+d_12*b)/(-f_1*real(x(i))*p_2/f_2+p_1-d_12);
        
        cos_theta = real(x(i));
        sin_theta = sqrt(1-real(x(i))^2);
        sin_alpha = sqrt(1/(cot_alpha^2+1));
        cos_alpha = sqrt(1-sin_alpha^2);
        
        if cot_alpha < 0
            cos_alpha = -cos_alpha;
        end

        C = [ d_12*cos_alpha*(sin_alpha*b+cos_alpha);
              cos_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha);
              sin_theta*d_12*sin_alpha*(sin_alpha*b+cos_alpha) ];
          
        C = P1 + N'*C;

        R = [ -cos_alpha -sin_alpha*cos_theta -sin_alpha*sin_theta;
              sin_alpha -cos_alpha*cos_theta -cos_alpha*sin_theta;
              0 -sin_theta cos_theta ];
        
        R = N'*R'*T;
        
        poses(1:3,(i-1)*4+1) = C;
        poses(1:3,(i-1)*4+2:(i-1)*4+4) = R;
    
    end
    
end
