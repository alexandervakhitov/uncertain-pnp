function [R,t,Xc,best_solution]=efficient_pnp_planar_gauss_lines(x3d_h,x2d_h,Xs, Xe, xs, xe, A, l2d)

% EFFICIENT_PNP Main Function to solve the PnP problem 
%       as described in:
%
%       Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua.
%       Accurate Non-Iterative O(n) Solution to the PnP Problem. 
%       In Proceedings of ICCV, 2007. 
%
%       x3d_h: homogeneous coordinates of the points in world reference
%       x2d_h: homogeneous position of the points in the image plane
%       A: intrincic camera parameters
%       R: Rotation of the camera system wrt world reference
%       T: Translation of the camera system wrt world reference
%       Xc: Position of the points in the camera reference
%
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


%=============================================================
% Modified by
% Xu Chi, Huazhong University of Science and Technology, China
%=============================================================


Xw = x3d_h(1:3, :);
U = x2d_h(1:2, :);

R= [];
t= [];

THRESHOLD_REPROJECTION_ERROR=5/800;%error in degrees of the basis formed by the control points. 
%If we have a larger error, we will compute the solution using a larger
%number of vectors in the kernel

%define control points in a world coordinate system (centered on the 3d
%points centroid)
Cw=[1 0 0;
    0 1 0;
    0 0 0];

%compute alphas (linear combination of the control points to represent the 3d
%points)
[Alph, Alph_c, ma] = epnp_orig.compute_alphas_planar(Xw, Cw);
[Alph_S, Alph_Sc, ms] = epnp_orig.compute_alphas_planar(Xs, Cw);

if (length(Xw) > 3)
    [Uac,Sac, ~] = svd(Alph_c*Alph_c');
    RA = Uac*Sac(1:3,1:3)*Uac'*Cw;
    ac_svals = sqrt(diag(Sac(1:3, 1:3))');
    SUac = diag(ac_svals)*Uac';
    temp_cw =  SUac*Cw;
    norm_cw = norm(temp_cw(:));
    Xw_mean = mean(Xw, 2);
else
    [Uac,Sac, ~] = svd(Alph_Sc*Alph_Sc');
    RA = Uac*Sac(1:3,1:3)*Uac'*Cw;
    ac_svals = sqrt(diag(Sac(1:3, 1:3))');
    SUac = diag(ac_svals)*Uac';
    temp_cw =  SUac*Cw;
    norm_cw = norm(temp_cw(:));
    Xs_mean = mean(Xs, 2);
end
Alph_E = epnp_orig.compute_alphas_planar(Xe, Cw);

%Compute M
% l2d = computeImageLines(xs, xe);
Km = epnp_lines.compute_M_ver3_planar(U, Alph,A, Alph_S, Alph_E, l2d);

%1.-Solve assuming dim(ker(M))=1. X=[Km_end];------------------------------
dim_kerM=1;
X1=Km(:,end);

if (length(Xw) <= 3)
    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X1, Alph_S, SUac, norm_cw);    
    [R,T] = epnp_lines.getrotTPlanar(Xs_mean, RA, Cc, ms);  %solve exterior orientation
    err(1) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);%reprojection_error_usingRT(Xw,U,R,T,A);
else
    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X1, Alph, SUac, norm_cw);
    [R,T] = epnp_lines.getrotTPlanar(Xw_mean, RA, Cc, ma);  %solve exterior orientation
    err(1) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);    
end

sol(1).Xc=Xc;
sol(1).Cc=Cc;
sol(1).R=R;
sol(1).T=T;
sol(1).error=err(1);
sol(1).sc = sc;
sol(1).betas = 1;


%2.-Solve assuming dim(ker(M))=2------------------------------------------
Km1=Km(:,end-1);
Km2=Km(:,end);

%control points distance constraint
D=compute_constraint_distance_2param(Km1,Km2);
dsq=define_distances_btw_control_points2();
betas_=inv(D'*D)*D'*dsq;
beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2))*sign(betas_(1));
X2=beta1*Km1+beta2*Km2;
if (norm(X2) == 0)
    X2 = ones(9,1);
end

if (length(Xw) <= 3)
    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X2, Alph_S, SUac, norm_cw);    
    [R,T] = epnp_lines.getrotTPlanar(Xs_mean, RA, Cc, ms);  %solve exterior orientation
    err(2) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);%reprojection_error_usingRT(Xw,U,R,T,A);
else
    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X2, Alph, SUac, norm_cw);
    [R,T] = epnp_lines.getrotTPlanar(Xw_mean, RA, Cc, ma);  %solve exterior orientation
    err(2) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);    
end

sol(2).Xc=Xc;
sol(2).Cc=Cc;
sol(2).R=R;
sol(2).T=T;
sol(2).error=err(2);
sol(2).sc = sc;
sol(2).betas = [beta1, beta2];

%3.-Solve assuming dim(ker(M))=3------------------------------------------
if min(err)>THRESHOLD_REPROJECTION_ERROR %just compute if we do not have good solution in the previus cases

    Km1=Km(:,end-2);
    Km2=Km(:,end-1);
    Km3=Km(:,end);

    %control points distance constraint
    D=compute_constraint_distance_3param(Km1,Km2,Km3);
    dsq=define_distances_btw_control_points2();
    D_=[D,-dsq];
    Kd=null(D_);
    for i= 1:4
        Kd(:,i)= Kd(:,i)/Kd(7,i);
    end
    Kd= Kd(1:6,:);
    P=compute_permutation_constraint3(Kd);

    lambdas_ = epnp_lines.kernel_noise(P,1);
    lambda(1)=sqrt(abs(lambdas_(1)));
    lambda(2)=sqrt(abs(lambdas_(5)))*sign(lambdas_(2));
    lambda(3)=sqrt(abs(lambdas_(8)))*sign(lambdas_(3));
    lambda(4)=sqrt(abs(lambdas_(10)))*sign(lambdas_(4));
    
    betas_=lambda(1)*Kd(:,1)+lambda(2)*Kd(:,2)+lambda(3)*Kd(:,3)+lambda(4)*Kd(:,4);
    beta1=sqrt(abs(betas_(1)));
    beta2=sqrt(abs(betas_(4)))*sign(betas_(2))*sign(betas_(1));
    beta3=sqrt(abs(betas_(6)))*sign(betas_(3))*sign(betas_(1));

    X3=beta1*Km1+beta2*Km2+beta3*Km3;

    if (length(Xw) <= 3)
        [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X3, Alph_S, SUac, norm_cw);    
        [R,T] = epnp_lines.getrotTPlanar(Xs_mean, RA, Cc, ms);  %solve exterior orientation
        err(3) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);%reprojection_error_usingRT(Xw,U,R,T,A);
    else
        [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor2(X3, Alph, SUac, norm_cw);
        [R,T] = epnp_lines.getrotTPlanar(Xw_mean, RA, Cc, ma);  %solve exterior orientation
        err(3) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);    
    end

    sol(3).Xc=Xc;
    sol(3).Cc=Cc;
    sol(3).R=R;
    sol(3).T=T;
    sol(3).error=err(3);
    sol(3).sc = sc;
    sol(3).betas = [beta1 beta2 beta3];

end

%4.-Gauss Newton Optimization------------------------------------------------------ 
[min_err,best_solution]=min(err);
Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
t=sol(best_solution).T;
Betas=sol(best_solution).betas;
sc=sol(best_solution).sc;
% Kernel=sol(best_solution).Kernel;
 
if best_solution==1
    Betas=[0,0,Betas];
elseif best_solution==2
    Betas=[0,Betas];
elseif best_solution==3
    Betas=[Betas];
end


Km1=Km(:,end-2);
Km2=Km(:,end-1);
Km3=Km(:,end);
Kernel=[Km1,Km2,Km3];


%refine the solution iterating over the betas
Beta0=Betas/sc;
if (length(Xw) <= 3)
    [Xc_opt,R_opt,T_opt,err_opt,iter] = epnp_lines.optimize_betas_planar_gauss_newton_lines(Kernel,...
        Cw,Beta0,Xs,Xe,xs,xe,A,Xw,U,Alph,Alph_S,Xs_mean, RA, ms);
else
    [Xc_opt,R_opt,T_opt,err_opt,iter] = epnp_lines.optimize_betas_planar_gauss_newton_lines(Kernel,...
        Cw,Beta0,Xs,Xe,xs,xe,A,Xw,U,Alph,Alph_S,Xw_mean, RA, ma);
    %(Km,...
    %Cw,Beta0,Xs,Xe,xs,xe,A,Xw, U,Alpha, Alpha_S, ptMean, RA, ma)
end

% err_opt = 1;
%Just update R,T,Xc if Gauss Newton improves results (which is almost
%always)
if err_opt<min_err    
    R=R_opt;
    t=T_opt;
    Xc=Xc_opt;
end

opt.Beta0=Beta0;
opt.Kernel=Kernel;
opt.iter=iter;




return

%4.-Solve assuming dim(ker(M))=4------------------------------------------
if min(err)>THRESHOLD_REPROJECTION_ERROR %just compute if we do not have good solution in the previus cases
    Km1=Km(:,end-3);
    Km2=Km(:,end-2);
    Km3=Km(:,end-1);
    Km4=Km(:,end);


    D=compute_constraint_distance_orthog_4param_9eq_10unk(Km1,Km2,Km3,Km4);
    dsq=define_distances_btw_control_points();
    lastcolumn=[-dsq',0,0,0]';
    D_=[D,lastcolumn];
    Kd=null(D_);

    P=compute_permutation_constraint4(Kd);
    lambdas_=kernel_noise(P,1);
    lambda(1)=sqrt(abs(lambdas_(1)));
    lambda(2)=sqrt(abs(lambdas_(6)))*sign(lambdas_(2))*sign(lambdas_(1));
    lambda(3)=sqrt(abs(lambdas_(10)))*sign(lambdas_(3))*sign(lambdas_(1));
    lambda(4)=sqrt(abs(lambdas_(13)))*sign(lambdas_(4))*sign(lambdas_(1));
    lambda(5)=sqrt(abs(lambdas_(15)))*sign(lambdas_(5))*sign(lambdas_(1));

    betass_=lambda(1)*Kd(:,1)+lambda(2)*Kd(:,2)+lambda(3)*Kd(:,3)+lambda(4)*Kd(:,4)+lambda(5)*Kd(:,5);
    beta1=sqrt(abs(betass_(1)));
    beta2=sqrt(abs(betass_(5)))*sign(betass_(2));
    beta3=sqrt(abs(betass_(8)))*sign(betass_(3));
    beta4=sqrt(abs(betass_(10)))*sign(betass_(4));
    X4=beta1*Km1+beta2*Km2+beta3*Km3+beta4*Km4;

    [Cc,Xc]=compute_norm_sign_scaling_factor(X4,Cw,Alph,Xw);
    
    [R,T]=getrotTPlanarPlanar(Xw,Xc);  %solve exterior orientation
    err(4)=reprojection_error_usingRT(Xw,U,R,T,A);

    sol(4).Xc=Xc;
    sol(4).Cc=Cc;
    sol(4).R=R;
    sol(4).T=T;
    sol(4).error=err(4);

end

[min_err,best_solution]=min(err);
Xc=sol(best_solution).Xc;
R=sol(best_solution).R;
t=sol(best_solution).T;

return

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

return

function M=compute_M2(U,Alph,A)

n=size(Alph,1); %number of 3d points

fu=A(1,1);
fv=A(2,2);
u0=A(1,3);
v0=A(2,3);

nrows_M=2*n;
ncols_M=9;
M=zeros(nrows_M,ncols_M);

for i=1:n
    a1=Alph(i,1);
    a2=Alph(i,2);
    a3=Alph(i,3);
    
    ui=U(i,1);
    vi=U(i,2);
    
    %generate submatrix M
    M_=[a1*fu, 0, a1*(u0-ui), a2*fu, 0, a2*(u0-ui), a3*fu, 0, a3*(u0-ui);
        0, a1*fv, a1*(v0-vi), 0, a2*fv, a2*(v0-vi), 0, a3*fv, a3*(v0-vi)];
    
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=M_;
     
end

return


function P=compute_constraint_distance_2param(m1,m2)

%redefine variables name, for compatibility with maple
m1_1=m1(1); 
m1_2=m1(2); 
m1_3=m1(3); 
m1_4=m1(4); 
m1_5=m1(5); 
m1_6=m1(6);
m1_7=m1(7); 
m1_8=m1(8); 
m1_9=m1(9); 

m2_1=m2(1); 
m2_2=m2(2); 
m2_3=m2(3); 
m2_4=m2(4); 
m2_5=m2(5); 
m2_6=m2(6);
m2_7=m2(7); 
m2_8=m2(8); 
m2_9=m2(9); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t7 = (m1_6 ^ 2);
t8 = (m1_4 ^ 2);
t9 = (m1_1 ^ 2);
t10 = (m1_5 ^ 2);
t11 = (m1_2 ^ 2);
t12 = (m1_3 ^ 2);
t17 = m1_4 * m2_4;
t18 = m1_1 * m2_1;
t19 = m1_5 * m2_5;
t22 = m1_2 * m2_2;
t23 = m1_6 * m2_6;
t25 = m1_3 * m2_3;
t26 = (-m2_6 * m1_3 - m1_4 * m2_1 - m2_4 * m1_1 + t17 + t18 + t19 - m1_5 * m2_2 - m2_5 * m1_2 + t22 + t23 - m1_6 * m2_3 + t25);
t29 = (m2_3 ^ 2);
t34 = (m2_4 ^ 2);
t35 = (m2_1 ^ 2);
t36 = (m2_5 ^ 2);
t37 = (m2_2 ^ 2);
t38 = (m2_6 ^ 2);
t44 = (m1_7 ^ 2);
t45 = (m1_8 ^ 2);
t46 = (m1_9 ^ 2);
t55 = m1_8 * m2_8;
t56 = m1_9 * m2_9;
t58 = m1_7 * m2_7;
t59 = (-m1_9 * m2_3 - m2_8 * m1_2 - m2_9 * m1_3 - m1_7 * m2_1 - m2_7 * m1_1 + t55 + t22 + t56 + t18 - m1_8 * m2_2 + t25 + t58);
t64 = (m2_8 ^ 2);
t65 = (m2_9 ^ 2);
t68 = (m2_7 ^ 2);

t113 = (-m1_9 * m2_6 - m2_9 * m1_6 + t55 + t23 + t17 + t56 + t58 - m1_7 * m2_4 - m2_7 * m1_4 - m1_8 * m2_5 - m2_8 * m1_5 + t19);
P(1,1) = -2 * m1_4 * m1_1 - 2 * m1_5 * m1_2 - 2 * m1_6 * m1_3 + t7 + t8 + t9 + t10 + t11 + t12;
P(1,2) = 2 * t26;
P(1,3) = -2 * m2_6 * m2_3 + t29 - 2 * m2_4 * m2_1 - 2 * m2_5 * m2_2 + t34 + t35 + t36 + t37 + t38;
P(2,1) = -2 * m1_7 * m1_1 + t12 - 2 * m1_9 * m1_3 + t44 + t45 + t46 - 2 * m1_8 * m1_2 + t9 + t11;
P(2,2) = 2 * t59;
P(2,3) = -2 * m2_8 * m2_2 - 2 * m2_9 * m2_3 + t64 + t65 - 2 * m2_7 * m2_1 + t29 + t68 + t37 + t35;
P(3,1) = -2 * m1_9 * m1_6 + t8 + t10 + t7 - 2 * m1_7 * m1_4 + t44 + t45 + t46 - 2 * m1_8 * m1_5;
P(3,2) = 2 * t113;
P(3,3) = -2 * m2_9 * m2_6 + t68 + t64 - 2 * m2_7 * m2_4 - 2 * m2_8 * m2_5 + t34 + t36 + t38 + t65;

return

function P=compute_constraint_distance_3param(m1,m2,m3)

%redefine variables name, for compatibility with maple
m1_1=m1(1); 
m1_2=m1(2); 
m1_3=m1(3); 
m1_4=m1(4); 
m1_5=m1(5); 
m1_6=m1(6);
m1_10=m1(7); 
m1_11=m1(8); 
m1_12=m1(9);

m2_1=m2(1); 
m2_2=m2(2); 
m2_3=m2(3); 
m2_4=m2(4); 
m2_5=m2(5); 
m2_6=m2(6);
m2_10=m2(7); 
m2_11=m2(8); 
m2_12=m2(9);

m3_1=m3(1); 
m3_2=m3(2); 
m3_3=m3(3); 
m3_4=m3(4); 
m3_5=m3(5); 
m3_6=m3(6);
m3_10=m3(7); 
m3_11=m3(8); 
m3_12=m3(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1 = (m1_2 ^ 2);
t4 = (m1_6 ^ 2);
t5 = (m1_3 ^ 2);
t6 = (m1_5 ^ 2);
t11 = (m1_4 ^ 2);
t12 = (m1_1 ^ 2);
t20 = m1_4 * m2_4;
t21 = m1_3 * m2_3;
t22 = m1_5 * m2_5;
t23 = m1_2 * m2_2;
t24 = m1_6 * m2_6;
t25 = m1_1 * m2_1;
t26 = (-m2_4 * m1_1 - m2_5 * m1_2 - m2_6 * m1_3 - m1_6 * m2_3 - m1_4 * m2_1 - m1_5 * m2_2 + t20 + t21 + t22 + t23 + t24 + t25);
t27 = m1_6 * m3_6;
t29 = m1_5 * m3_5;
t30 = m1_4 * m3_4;
t33 = m1_3 * m3_3;
t35 = m1_1 * m3_1;
t38 = m1_2 * m3_2;
t39 = (t27 - m1_6 * m3_3 + t29 + t30 - m1_4 * m3_1 - m3_6 * m1_3 + t33 - m3_5 * m1_2 + t35 - m3_4 * m1_1 - m1_5 * m3_2 + t38);
t40 = (m2_4 ^ 2);
t41 = (m2_2 ^ 2);
t42 = (m2_5 ^ 2);
t43 = (m2_1 ^ 2);
t44 = (m2_6 ^ 2);
t45 = (m2_3 ^ 2);
t53 = m2_4 * m3_4;
t56 = m2_5 * m3_5;
t57 = m2_2 * m3_2;
t60 = m2_1 * m3_1;
t62 = m2_6 * m3_6;
t63 = m2_3 * m3_3;
t65 = (t53 - m2_4 * m3_1 - m3_4 * m2_1 + t56 + t57 - m2_5 * m3_2 - m2_6 * m3_3 + t60 - m3_6 * m2_3 + t62 + t63 - m3_5 * m2_2);
t66 = (m3_5 ^ 2);
t69 = (m3_4 ^ 2);
t70 = (m3_3 ^ 2);
t71 = (m3_2 ^ 2);
t72 = (m3_6 ^ 2);
t75 = (m3_1 ^ 2);
t141 = (m1_10 ^ 2);
t142 = (m1_11 ^ 2);
t143 = (m1_12 ^ 2);
t151 = m1_10 * m2_10;
t152 = m1_12 * m2_12;
t154 = m1_11 * m2_11;
t158 = (-m2_10 * m1_1 - m2_12 * m1_3 + t151 + t23 + t25 + t152 + t21 - m1_10 * m2_1 + t154 - m1_11 * m2_2 - m2_11 * m1_2 - m1_12 * m2_3);
t160 = m1_12 * m3_12;
t164 = m1_10 * m3_10;
t165 = m1_11 * m3_11;
t168 = (-m3_10 * m1_1 + t160 - m3_11 * m1_2 + t38 + t33 - m1_12 * m3_3 - m3_12 * m1_3 + t164 + t35 + t165 - m1_11 * m3_2 - m1_10 * m3_1);
t169 = (m2_12 ^ 2);
t170 = (m2_10 ^ 2);
t171 = (m2_11 ^ 2);
t179 = m2_10 * m3_10;
t181 = m2_12 * m3_12;
t185 = m2_11 * m3_11;
t188 = (t57 + t60 + t179 - m2_10 * m3_1 + t181 - m2_12 * m3_3 - m3_12 * m2_3 - m3_10 * m2_1 + t63 + t185 - m2_11 * m3_2 - m3_11 * m2_2);
t191 = (m3_12 ^ 2);
t192 = (m3_10 ^ 2);
t193 = (m3_11 ^ 2);
t254 = (-m2_11 * m1_5 + t154 + t151 - m1_12 * m2_6 + t20 - m1_10 * m2_4 - m2_12 * m1_6 - m2_10 * m1_4 - m1_11 * m2_5 + t152 + t24 + t22);
t261 = (t30 - m3_12 * m1_6 - m1_10 * m3_4 - m1_11 * m3_5 + t160 + t27 + t164 + t165 - m3_11 * m1_5 - m3_10 * m1_4 - m1_12 * m3_6 + t29);
t275 = (-m3_10 * m2_4 + t56 - m2_10 * m3_4 + t62 - m3_12 * m2_6 - m2_11 * m3_5 + t53 - m3_11 * m2_5 - m2_12 * m3_6 + t179 + t181 + t185);
% 12
P(1,1) = t1 - 2 * m1_4 * m1_1 + t4 + t5 + t6 - 2 * m1_5 * m1_2 - 2 * m1_6 * m1_3 + t11 + t12;
P(1,2) = 2 * t26;
P(1,3) = 2 * t39;
P(1,4) = t40 + t41 + t42 + t43 + t44 + t45 - 2 * m2_5 * m2_2 - 2 * m2_6 * m2_3 - 2 * m2_4 * m2_1;
P(1,5) = 2 * t65;
P(1,6) = t66 - 2 * m3_4 * m3_1 + t69 + t70 + t71 + t72 - 2 * m3_6 * m3_3 + t75 - 2 * m3_5 * m3_2;
% 14
P(2,1) = -2 * m1_11 * m1_2 + t141 + t142 + t12 + t1 + t5 + t143 - 2 * m1_10 * m1_1 - 2 * m1_12 * m1_3;
P(2,2) = 2 * t158;
P(2,3) = 2 * t168;
P(2,4) = t169 + t41 + t43 + t45 + t170 + t171 - 2 * m2_12 * m2_3 - 2 * m2_10 * m2_1 - 2 * m2_11 * m2_2;
P(2,5) = 2 * t188;
P(2,6) = t71 - 2 * m3_12 * m3_3 + t75 + t191 + t70 + t192 + t193 - 2 * m3_10 * m3_1 - 2 * m3_11 * m3_2;
% 24
P(3,1) = t4 + t143 + t11 - 2 * m1_12 * m1_6 - 2 * m1_11 * m1_5 - 2 * m1_10 * m1_4 + t6 + t141 + t142;
P(3,2) = 2 * t254;
P(3,3) = 2 * t261;
P(3,4) = t170 + t171 + t169 - 2 * m2_10 * m2_4 - 2 * m2_11 * m2_5 - 2 * m2_12 * m2_6 + t40 + t42 + t44;
P(3,5) = 2 * t275;
P(3,6) = t69 + t66 + t72 - 2 * m3_12 * m3_6 - 2 * m3_10 * m3_4 + t193 - 2 * m3_11 * m3_5 + t192 + t191;

return

function dsq=define_distances_btw_control_points2()

%relative coordinates of the control points
c1=[1,0,0];
c2=[0,1,0];
c3=[0,0,0];

d12=(c1(1)-c2(1))^2 + (c1(2)-c2(2))^2 + (c1(3)-c2(3))^2;
d13=(c1(1)-c3(1))^2 + (c1(2)-c3(2))^2 + (c1(3)-c3(3))^2;
d23=(c2(1)-c3(1))^2 + (c2(2)-c3(2))^2 + (c2(3)-c3(3))^2;

dsq=[d12,d13,d23]';

return

function K=compute_permutation_constraint3(V)

%[B11,B12,...,B33]=lambda1*v1+lambda2*v2+lambda3*v3

N=size(V,2); %dimension of the kernel
n=3; %dimension of Bij
idx=[1 2 3; 2 4 5; 3 5 6];

%1.-Generation of the first set of equations Bii.Bjj=Bij.Bii  (n(n-1)/2 eqs).
nrowsK=n*(n-1)/2+n*(n-1)*n/2;
ncolsK=N*(N+1)/2;
K=zeros(nrowsK,ncolsK);

t=1;
for i=1:n
    for j=i+1:n
        offset=1;
        for a=1:N
            for b=a:N
                if a==b
                    K(t,offset)=V(idx(i,i),a)*V(idx(j,j),a)-V(idx(i,j),a)*V(idx(i,j),a);
                else
                    K(t,offset)=V(idx(i,i),a)*V(idx(j,j),b)-V(idx(i,j),a)*V(idx(i,j),b)+...
                                V(idx(i,i),b)*V(idx(j,j),a)-V(idx(i,j),b)*V(idx(i,j),a);
                end
                offset=offset+1;
            end
            
        end
        t=t+1;
        %fprintf('t:%d\t offset:%d\n',t,offset);
    end
end


for k=1:n
    for j=k:n
        for i=1:n
            if (i~=j & i~=k)
                offset=1;
                for a=1:N
                    for b=a:N
                        if a==b
                            K(t,offset)=V(idx(i,j),a)*V(idx(i,k),a)-V(idx(i,i),a)*V(idx(j,k),a);
                        else
                            K(t,offset)=V(idx(i,j),a)*V(idx(i,k),b)-V(idx(i,i),a)*V(idx(j,k),b)+...
                                        V(idx(i,j),b)*V(idx(i,k),a)-V(idx(i,i),b)*V(idx(j,k),a);
                        end
                        offset=offset+1;
                    end
                    
                end
                t=t+1;
                %fprintf('t:%d\t offset:%d\n',t,offset);
            end
        end
    end
end
                
return
