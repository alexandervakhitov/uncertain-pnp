function [R, T]=getrotT(wcent, RA, Cc, ma)
  
% This routine solves the exterior orientation problem for a point cloud
%  given in both camera and world coordinates. 
  
% wpts = 3D points in arbitrary reference frame
% cpts = 3D points in camera reference frame

M = Cc'*RA;
[U S V]=svd(M);
R=U*V';
if det(R)<0
  R=-R;
end
T = Cc'*ma-R*wcent;
% 
end
