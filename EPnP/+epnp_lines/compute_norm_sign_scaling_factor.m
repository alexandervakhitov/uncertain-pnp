function [Cc,Xc,sc]=compute_norm_sign_scaling_factor(X1,Alph,SUac, norm_cw)
 
% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 


n=size(Alph, 2); %number of data points

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(4,3);
for i=1:4
    Cc_(i,:)=X1(3*i-2:3*i);
end

% %position of reference points in camera coordinates
% Xc_=Cc_'*Alph;
% 
% %compute distances in world coordinates w.r.t. the centroid
% centr_w = mean(Xw, 2);
% centroid_w = repmat(centr_w,[1,n]);
% tmp1 = Xw-centroid_w;
% dist_w = sqrt(sum(tmp1.^2,2));
% 
% %compute distances in camera coordinates w.r.t. the centroid
% centr_c = mean(Xc_, 2);
% centroid_c = repmat(centr_c,[1, n]);
% tmp2 = Xc_-centroid_c;
% dist_c = sqrt(sum(tmp2.^2,2));
%  
% %least squares solution for the scale factor
% sc=1/(inv(dist_c'*dist_c)*dist_c'*dist_w);
% 
% %scale position of the control points
% Cc=Cc_/sc;

%rescaled position of the reference points



temp_cc =  SUac*Cc_;
norm_cc = norm(temp_cc(:));
sc = norm_cc/norm_cw;
Cc = 1/sc*Cc_;
Xc = Cc'*Alph;

%change the sign if necessary. z negative is no possible in camera
%coordinates
neg_z=find(Xc(3, :)<0);
if size(neg_z,1)>=1
    sc=-sc;
    Xc=Xc*(-1);
end




        
        
