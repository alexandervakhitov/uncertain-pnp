function [Cc,Xc,sc]=compute_norm_sign_scaling_factor2(X1,Alph,SUac, norm_cw)
 

%Km will be a scaled solution. In order to find the scale parameter we
%impose distance constraints between the reference points

%scaled position of the control points in camera coordinates
Cc_=zeros(3,3);
for i=1:3
    Cc_(i,:)=X1(3*i-2:3*i);
end

temp_cc =  SUac*Cc_;
norm_cc = norm(temp_cc(:));
sc = norm_cc/norm_cw;
Cc = 1/sc*Cc_;
Xc = Cc'*Alph;

% change the sign if necessary. z negative is no possible in camera
% coordinates
neg_z=find(Xc(3,:)<0);
if size(neg_z,2)>=1
    sc=-sc;
    Xc=Xc*(-1);
    Cc=Cc*(-1);
end



return
