function [Cc,Xctot,sc]=compute_norm_sign_scaling_factor_l_pt(X1,Alph_p, Alph_tot, norm_cw, SUac_tot)
 
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


%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(4,3);
for i=1:4
    Cc_(i,:)=X1(3*i-2:3*i);
end


temp_cc_tot =  SUac_tot*Cc_;
% temp_cc_l =  SUac_l*Cc_;
% temp_cc = [temp_cc_p; temp_cc_l];
norm_cc = norm(temp_cc_tot(:));
sc = norm_cc/norm_cw;
Cc = 1/sc*Cc_;
Xc = Cc'*Alph_p;
Xctot = Cc'*Alph_tot;
%change the sign if necessary. z negative is no possible in camera
%coordinates
if (size(Xc, 1) > 0)
    neg_z=find(Xc(3, :)<0);
    if size(neg_z,1)>=1
        sc=-sc;
        Xc=Xc*(-1);    
        Xctot = -Xctot;
    end
end




        
        
