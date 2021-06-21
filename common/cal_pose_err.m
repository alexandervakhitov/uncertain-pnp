function [y  y1]= cal_pose_err(T1, T2)

R1= T1(1:3,1:3);
[U,S,V] = svd(R1);
R1 = U*V';
R2= T2(1:3,1:3);

% X1= R1(:,1); X2= R2(:,1);
% Y1= R1(:,2); Y2= R2(:,2);
% Z1= R1(:,3); Z2= R2(:,3);
% 
% exyz= [X1'*X2 Y1'*Y2 Z1'*Z2];
% exyz(exyz>1)= 1;
% exyz(exyz<-1)= -1;
% 
% y(1)= max(abs(acos(exyz)))*180/pi;
dR = R1'*R2;
y(1) = acos(0.5*(trace(dR)-1))*180/pi;

% q1 = Matrix2Quaternion(R1);
% q2 = Matrix2Quaternion(R2);
% 
% y1(1) = norm(q1-q2)/norm(q2)*100;
% 
% if isnan(y(1))
%     txt;
% end

y(2)= norm(T1(1:3,4)-T2(1:3,4))/(norm(T2(1:3,4)) + 1e-3)*100;
y1(2) = y(2);
y= abs(y);

y(3) = norm(-R1' * T1(1:3, 4) + R2' * T2(1:3, 4));

end


