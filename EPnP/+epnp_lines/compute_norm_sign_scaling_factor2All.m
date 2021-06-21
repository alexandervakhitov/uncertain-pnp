function [Cc,Xctot,sc]=compute_norm_sign_scaling_factor2All(X1,Alph_p, Alph_tot, norm_cw, SUac_tot)
 

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(3,3);
for i=1:3
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
if (size(Xc, 2) > 0)
    neg_z = find(Xc(3, :)<0);
else
    neg_z = find(Xctot(3, :) < 0);
end

if size(neg_z,2)>=1
    sc=-sc;
    Xc=Xc*(-1);    
    Xctot = -Xctot;
    Cc = -Cc;
end


return
