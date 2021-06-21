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

function R2 = Rodrigues2(R1)

[r,c] = size(R1);

%% Rodrigues Rotation Vector to Rotation Matrix
if ((r == 3) && (c == 1)) || ((r == 1) && (c == 3))
    wx = [  0   -R1(3)  R1(2);
           R1(3)   0   -R1(1);
          -R1(2)  R1(1)   0   ];
      
    omega_norm = sqrt(R1(1)^2 + R1(2)^2 + R1(3)^2);
    
    if (omega_norm < eps)
        R2 = eye(3);
    else
        R2 = eye(3) + ...
             sin(omega_norm)/omega_norm*wx + ...
             (1-cos(omega_norm))/omega_norm^2*wx^2;
    end

%% Rotation Matrix to Rodrigues Rotation Vector
elseif (r == 3) && (c == 3)
    w_norm = acos((trace(R1)-1)/2);
    if (w_norm < eps)
        R2 = [0 0 0]';
    else
        R2 = 1/(2*sin(w_norm)) * ...
            [R1(3,2)-R1(2,3);R1(1,3)-R1(3,1);R1(2,1)-R1(1,2)]*w_norm;
    end
end


