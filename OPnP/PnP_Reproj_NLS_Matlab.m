function [Rr tr] = PnP_Reproj_NLS_Matlab(U,u,R0,t0)
%to refine the column-triplet by using nonlinear least square
%it's much better than the fminunc used by CVPR12.
%there are three fractional formulations, anyone is equivalently good. 
%there are two constrained formulations, yet neither is good. 
%the reason is that: eliminating the unkown scale. ---> similar to
%reprojection error and homogeneous error

%transform matrix to unit quaternion
q = matrix2quaternion(R0);

%nonlinear least square parameter setting
options = optimset('Algorithm','trust-region-reflective','Jacobian','on','DerivativeCheck','off','Display','off');

%call matlab lsqnonlin
[var,resnorm] = lsqnonlin(@(var)valuegradient(var,U,u),[q(:);t0(:)],[],[],options);

%denormalize data
qr = var(1:4); tr = var(5:end);
normqr = norm(qr);
tr = tr/(normqr)^2;
qr = qr/normqr;

%trasform back into matrix form
Rr = quaternion2matrix(qr);
end

%formulation 1
function [fval grad] = valuegradient(var,U,u)
npt = size(U,2);
fval = zeros(2*npt,1);
grad = zeros(2*npt,7);

a = var(1); b = var(2); c = var(3); d = var(4);
R = [a^2+b^2-c^2-d^2     2*b*c-2*a*d     2*b*d+2*a*c
     2*b*c+2*a*d         a^2-b^2+c^2-d^2 2*c*d-2*a*b
     2*b*d-2*a*c         2*c*d+2*a*b   a^2-b^2-c^2+d^2];
t1 = var(5); t2 = var(6); t3 = var(7);
for i = 1:npt
    vec = R*U(:,i) + [t1;t2;t3];
    fval(2*(i-1)+1) = u(1,i) - vec(1)/vec(3);
    fval(2*i) =  u(2,i) - vec(2)/vec(3);
    
    temp1 = [2*[a -d c;b c d; -c b a; -d -a b]*U(:,i); 1; 0; 0].';
    temp2 = [2*[d a -b;c -b -a;b c d; a -d c]*U(:,i); 0; 1; 0].';
    temp3 = [2*[-c b a;d a -b;-a d -c; b c d]*U(:,i); 0; 0; 1].';
    
    grad(2*(i-1)+1,:) = (-temp1*vec(3) + temp3*vec(1))/(vec(3)^2);
    grad(2*i,:) = (-temp2*vec(3) + temp3*vec(2))/(vec(3)^2);
end
end

