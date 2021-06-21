%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       A. Vakhitov, L. Ferraz, A. Agudo, F. Moreno-Noguer
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       CVPR 2021
%
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R,T,Xc,best_solution,opt]=efficient_pnpu_gauss_fast(x3d_h, x2d_h, A,...
                                                              sigmas3D, sigmas2D, ...
                                                              Xs, Xe, l, ...
                                                              Sigmas3DLines, ...
                                                              Sigmas2DLines, ...
                                                              Rest, test)

Xw = x3d_h(:,1:3);
U = x2d_h(:,1:2);

THRESHOLD_REPROJECTION_ERROR=0.1;%error in degrees of the basis formed by the control points. 
%If we have a larger error, we will compute the solution using a larger
%number of vectors in the kernel

%define control points in a world coordinate system (centered on the 3d
%points centroid)

Cw=define_control_points_pca_unc_fast(Xw', sigmas3D);    

%compute alphas (linear combination of the control points to represent the 3d
%points)
Alph=compute_alphas(Xw,Cw);
[Alph_S, ~, ~] = epnp_lines.compute_alphas(Xs, Cw);
[Alph_E, ~, ~] = epnp_lines.compute_alphas(Xe, Cw);
Alph_S = Alph_S';
Alph_E = Alph_E';
%Compute M
M = epnp_unc.compute_M_uncertain_lines(U,Alph,A,sigmas3D,sigmas2D, Cw, Rest, ...
    test, Alph_S, Alph_E, l, Sigmas3DLines, Sigmas2DLines);
%Compute kernel M
Km = kernel_noise(M,4); %in matlab we have directly the funcion km=null(M);

%1.-Solve assuming dim(ker(M))=1. X=[Km_end];------------------------------
dim_kerM=1;
X1=Km(:,end);
[Cc,Xc,sc]=compute_norm_sign_scaling_factor_ext(X1,Cw,Alph,Xw);

[R,T]=getrotT(Xw,Xc);  %solve exterior orientation
err(1)=reprojection_error_usingRT(Xw,U,R,T,A);

sol(1).Xc=Xc;
sol(1).Cc=Cc;
sol(1).R=R;
sol(1).T=T;
sol(1).error=err(1);
sol(1).betas=[1];
sol(1).sc=sc;
sol(1).Kernel=X1;


%2.-Solve assuming dim(ker(M))=2------------------------------------------
Km1=Km(:,end-1);
Km2=Km(:,end);

%control points distance constraint
D=compute_constraint_distance_2param_6eq_3unk(Km1,Km2);
if true
    dsq=define_distances_btw_control_points_pca(Cw);
else
    dsq=define_distances_btw_control_points();
end
betas_=(D'*D)\D'*dsq;
beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2))*sign(betas_(1));
X2=beta1*Km1+beta2*Km2;

[Cc,Xc,sc]=compute_norm_sign_scaling_factor_ext(X2,Cw,Alph,Xw);

[R,T]=getrotT(Xw,Xc);  %solve exterior orientation
err(2)=reprojection_error_usingRT(Xw,U,R,T,A);

sol(2).Xc=Xc;
sol(2).Cc=Cc;
sol(2).R=R;
sol(2).T=T;
sol(2).error=err(2);
sol(2).betas=[beta1,beta2];
sol(2).sc=sc;
sol(2).Kernel=[Km1,Km2];

%3.-Solve assuming dim(ker(M))=3------------------------------------------
if min(err)>THRESHOLD_REPROJECTION_ERROR %just compute if we do not have good solution in the previus cases

    Km1=Km(:,end-2);
    Km2=Km(:,end-1);
    Km3=Km(:,end);

    %control points distance constraint
    D=compute_constraint_distance_3param_6eq_6unk(Km1,Km2,Km3);
    if true
        dsq=define_distances_btw_control_points_pca(Cw);
    else
        dsq=define_distances_btw_control_points();
    end

    betas_=inv(D)*dsq;
    beta1=sqrt(abs(betas_(1)));
    beta2=sqrt(abs(betas_(4)))*sign(betas_(2))*sign(betas_(1));
    beta3=sqrt(abs(betas_(6)))*sign(betas_(3))*sign(betas_(1));

    X3=beta1*Km1+beta2*Km2+beta3*Km3;

    [Cc,Xc,sc]=compute_norm_sign_scaling_factor_ext(X3,Cw,Alph,Xw);
  
    [R,T]=getrotT(Xw,Xc);  %solve exterior orientation
    err(3)=reprojection_error_usingRT(Xw,U,R,T,A);

    sol(3).Xc=Xc;
    sol(3).Cc=Cc;
    sol(3).R=R;
    sol(3).T=T;
    sol(3).error=err(3);
    sol(3).betas=[beta1,beta2,beta3];
    sol(3).sc=sc;
    sol(3).Kernel=[Km1,Km2,Km3];

end



%4.-Solve assuming dim(ker(M))=4------------------------------------------
if (false)
%if min(err)>THRESHOLD_REPROJECTION_ERROR %just compute if we do not have good solution in the previus cases
    Km1=Km(:,end-3);
    Km2=Km(:,end-2);
    Km3=Km(:,end-1);
    Km4=Km(:,end);


    D=compute_constraint_distance_orthog_4param_9eq_10unk(Km1,Km2,Km3,Km4);
        
    if true
        dsq=define_distances_btw_control_points_pca(Cw);
    else
        dsq=define_distances_btw_control_points();
    end

    lastcolumn=[-dsq',0,0,0]';
    D_=[D,lastcolumn];
    Kd=null(D_);

    P=compute_permutation_constraint4(Kd);        
    lambdas_=kernel_noise(P,1);            
    lambda(1)=sqrt(abs(lambdas_(1)));
    lambda(2)=sqrt(abs(lambdas_(6)))*sign(lambdas_(2))*sign(lambdas_(1));
    lambda(3)=sqrt(abs(lambdas_(10)))*sign(lambdas_(3))*sign(lambdas_(1));
%     size(lambdas_)
    lambda(4)=sqrt(abs(lambdas_(13)))*sign(lambdas_(4))*sign(lambdas_(1));
    lambda(5)=sqrt(abs(lambdas_(15)))*sign(lambdas_(5))*sign(lambdas_(1));

    betass_=lambda(1)*Kd(:,1)+lambda(2)*Kd(:,2)+lambda(3)*Kd(:,3)+lambda(4)*Kd(:,4)+lambda(5)*Kd(:,5);
    beta1=sqrt(abs(betass_(1)));
    beta2=sqrt(abs(betass_(5)))*sign(betass_(2));
    beta3=sqrt(abs(betass_(8)))*sign(betass_(3));
    beta4=sqrt(abs(betass_(10)))*sign(betass_(4));
    X4=beta1*Km1+beta2*Km2+beta3*Km3+beta4*Km4;

    [Cc,Xc,sc]=compute_norm_sign_scaling_factor_ext(X4,Cw,Alph,Xw);
    
    [R,T]=getrotT(Xw,Xc);  %solve exterior orientation
    err(4)=reprojection_error_usingRT(Xw,U,R,T,A);

    sol(4).Xc=Xc;
    sol(4).Cc=Cc;
    sol(4).R=R;
    sol(4).T=T;
    sol(4).error=err(4);
    sol(4).betas=[beta1,beta2,beta3,beta4];
    sol(4).sc=sc;
    sol(4).Kernel=[Km1,Km2,Km3,Km4];
end


%5.-Gauss Newton Optimization------------------------------------------------------ 
[min_err,best_solution]=min(err);
Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
T=sol(best_solution).T;
Betas=sol(best_solution).betas;
sc=sol(best_solution).sc;
Kernel=sol(best_solution).Kernel;
 
if best_solution==1
    Betas=[0,0,0,Betas];
elseif best_solution==2
    Betas=[0,0,Betas];
elseif best_solution==3
    Betas=[0,Betas];
end

Km1=Km(:,end-3);
Km2=Km(:,end-2);
Km3=Km(:,end-1);
Km4=Km(:,end);
Kernel=[Km1,Km2,Km3,Km4];


%refine the solution iterating over the betas
Beta0=Betas/sc;

tic;
[Xc_opt,R_opt,T_opt,err_opt,iter]=optimize_betas_gauss_newton(Kernel,Cw,Beta0,Alph,Xw,U,A);
tco = toc;
% fprintf('epnpu m %f / opt %f \n', tcm, tco);
%Just update R,T,Xc if Gauss Newton improves results (which is almost
%always)
if err_opt<min_err    
    R=R_opt;
    T=T_opt;
    Xc=Xc_opt;
end

opt.Beta0=Beta0;
opt.Kernel=Kernel;
opt.iter=iter;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R, T]=getrotT(wpts,cpts)
  
% This routine solves the exterior orientation problem for a point cloud
%  given in both camera and world coordinates. 
  
% wpts = 3D points in arbitrary reference frame
% cpts = 3D points in camera reference frame
  
n=size(wpts,1);
M=zeros(3);

ccent=mean(cpts);
wcent=mean(wpts);

for i=1:3
  cpts(:,i)=cpts(:,i)-ccent(i)*ones(n,1);
  wpts(:,i)=wpts(:,i)-wcent(i)*ones(n,1);
end
for i=1:n
   M=M+cpts(i,:)'*wpts(i,:);
end
[U S V]=svd(M);
R=U*V';
if det(R)<0
  R=-R;
end
T=ccent'-R*wcent';
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)

%clear all; close all; load reprojection_error_usingRT;
n=size(Xw,1);

P=A*[R,T];
Xw_h=[Xw,ones(n,1)];
Urep_=(P*Xw_h')';

%project reference points into the image plane
Urep=zeros(n,2);
Urep(:,1)=Urep_(:,1)./Urep_(:,3);
Urep(:,2)=Urep_(:,2)./Urep_(:,3);


%reprojection error
err_=sqrt((U(:,1)-Urep(:,1)).^2+(U(:,2)-Urep(:,2)).^2);
err=sum(err_)/n;