function [R,t] = gOp_interface(U,u)
%make sure that sedumi is in the path. 
if exist('sedumi') < 1
    disp('Please install SeDuMi first for the SDP solution.')
    return;
end

opt.acc = 1e-8; 
opt.methode ='choose best';

if size(u,1) == 2
    u = [u;ones(1,size(u,2))];
end

% call gOp
[R,t] = gOp(u,U,opt);

end