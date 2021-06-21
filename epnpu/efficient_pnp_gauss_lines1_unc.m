function [R,T,Xc,best_solution,opt,t1,t2,t3,t4]=efficient_pnp_gauss_lines1_unc(x3d_h,x2d_h,Xs, Xe, xs, xe, A, l2d, Sigmas3D, Sigmas2D, SigmasLines, SigmasLines2D, Rest, test)

% EFFICIENT_PNP_GAUSS Main Function to solve the PnP problem 
%       as described in:
%
%       Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua.
%       Accurate Non-Iterative O(n) Solution to the PnP Problem. 
%       In Proceedings of ICCV, 2007. 
%
%       Note: In this version of the software we perform a final
%       optimization using Gauss-Newton,which is not described in the
%       paper.
%
%       x3d_h: homogeneous coordinates of the points in world reference
%       x2d_h: homogeneous position of the points in the image plane
%       A: intrincic camera parameters
%       R: Rotation of the camera system wrt world reference
%       T: Translation of the camera system wrt world reference
%       Xc: Position of the points in the camera reference
%       best solution: dimension of the kernel for the best solution
%                     (before applying Gauss Newton).
%       opt: some parameters of the optimization process
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
% Francesc Moreno-Noguer, CVLab-EPFL, October 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

t0 = tic;

if (nargin < 9)
    lnc = 1;
end


Xw=x3d_h(1:3,:);
U=x2d_h(1:2, :);

THRESHOLD_REPROJECTION_ERROR=1;%error in degrees of the basis formed by the control points. 
%If we have a larger error, we will compute the solution using a larger
%number of vectors in the kernel

%define control points in a world coordinate system (centered on the 3d
%points centroid)
Cw=define_control_points_pca_unc(Xw, Sigmas3D);

%compute alphas (linear combination of the control points to represent the 3d
%points)

t0 = tic;

[Alph, Alph_c, ma] = epnp_lines.compute_alphas(Xw, Cw);
[Alph_S, Alph_Sc, ms] = epnp_lines.compute_alphas(Xs, Cw);

Xtot = [Xw Xs Xe];
[Alph_T, Alph_Tc, m_tot] = epnp_lines.compute_alphas(Xtot, Cw);

[Utotc,Stotc, ~] = svd(Alph_Tc*Alph_Tc');
RA_tot = Utotc*Stotc(1:4,1:4)*Utotc'*Cw;
atotc_svals = sqrt(diag(Stotc(1:4, 1:4))');
SUac_tot = diag(atotc_svals)*Utotc';
temp_cw_tot =  SUac_tot*Cw;
norm_cw_tot = norm(temp_cw_tot(:));
XT_mean = mean(Xtot, 2);


% [Uac,Sac, ~] = svd(Alph_c*Alph_c');
% RA_p = Uac*Sac(1:4,1:4)*Uac'*Cw;
% ac_svals = sqrt(diag(Sac(1:4, 1:4))');
% SUac_p = diag(ac_svals)*Uac';
% temp_cw_p =  SUac*Cw;
% % norm_cwp = norm(temp_cw(:));
% Xw_mean = mean(Xw, 2);
% 
% [Uac,Sac, ~] = svd(Alph_Sc*Alph_Sc');
% RA_l = Uac*Sac(1:4,1:4)*Uac'*Cw;
% ac_svals = sqrt(diag(Sac(1:4, 1:4))');
% SUac_l = diag(ac_svals)*Uac';
% temp_cw_l =  SUac*Cw;
% % norm_cwl = norm(temp_cw(:));
% Xs_mean = mean(Xs, 2);
% norm_cw = norm([temp_cw_p; temp_cw_l]);



Alph_E = epnp_lines.compute_alphas(Xe, Cw);


t4 = toc(t0);
t0 = tic;
%Compute M
%Km = compute_M_ver3_unc(U, Alph,A, Alph_S, Alph_E, l2d, Sigmas3D, Sigmas2D, SigmasLines, SigmasLines2D, Cw, Rest, test)
Km = compute_M_ver3_unc(U, Alph, A, Alph_S, Alph_E, l2d, Sigmas3D, Sigmas2D, SigmasLines, SigmasLines2D, Cw, Rest, test);
t1 = toc(t0);

t0 = tic;
%Compute kernel M
% Km = epnp_orig.kernel_noise(M,4); %in matlab we have directly the funcion km=null(M);

%1.-Solve assuming dim(ker(M))=1. X=[Km_end];------------------------------
dim_kerM=1;
X1=Km(:,end);
[Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor_l_pt(X1, Alph, Alph_T, norm_cw_tot, SUac_tot);
[R,T] = epnp_lines.getrotT(XT_mean, RA_tot, Cc, m_tot);
err(1) = epnp_lines.reprojection_error_usingRTAll(Xw, U, Xs,Xe,xs,xe,R,T,A);

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
D = epnp_orig.compute_constraint_distance_2param_6eq_3unk(Km1,Km2);
dsq = epnp_orig.define_distances_btw_control_points();
%AV:inv -> pinv
betas_=pinv(D'*D)*D'*dsq;
beta1=sqrt(abs(betas_(1)));
beta2=sqrt(abs(betas_(3)))*sign(betas_(2))*sign(betas_(1));
X2=beta1*Km1+beta2*Km2;
if (norm(X2) == 0)
    X2 = Km1;
end

% if (length(Xw) < 3)    
%     [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor(X2,Alph_S, SUac, norm_cw);
%     [R,T] = epnp_lines.getrotT(Xs_mean, RA, Cc, ms); %solve exterior orientation
%     err(2) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);%reprojection_error_usingRT(Xw,U,R,T,A);    
% else
%     [Cc,Xc,sc]= epnp_lines.compute_norm_sign_scaling_factor(X2,Alph, SUac, norm_cw);
%     [R,T] = epnp_lines.getrotT(Xw_mean, RA, Cc, ma);  %solve exterior orientation
%     err(2) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);    
% end

[Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor_l_pt(X2, Alph, Alph_T, norm_cw_tot, SUac_tot);
[R,T] = epnp_lines.getrotT(XT_mean, RA_tot, Cc, m_tot);
err(2) = epnp_lines.reprojection_error_usingRTAll(Xw, U, Xs,Xe,xs,xe,R,T,A);

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
%     fprintf('solve 3\n');
    Km1=Km(:,end-2);
    Km2=Km(:,end-1);
    Km3=Km(:,end);

    %control points distance constraint
    D=epnp_orig.compute_constraint_distance_3param_6eq_6unk(Km1,Km2,Km3);
    dsq=epnp_orig.define_distances_btw_control_points();
    betas_=inv(D)*dsq;
    beta1=sqrt(abs(betas_(1)));
    beta2=sqrt(abs(betas_(4)))*sign(betas_(2))*sign(betas_(1));
    beta3=sqrt(abs(betas_(6)))*sign(betas_(3))*sign(betas_(1));

    X3=beta1*Km1+beta2*Km2+beta3*Km3;
%     if (length(Xw) <= 3)
%         [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor(X3,Alph_S, SUac, norm_cw);
%         [R,T] = epnp_lines.getrotT(Xs_mean, RA, Cc, ms); %solve exterior orientation
%         err(3) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);
%     else
%         [Cc,Xc,sc]=epnp_lines.compute_norm_sign_scaling_factor(X3,Alph, SUac, norm_cw);
%         [R,T] = epnp_lines.getrotT(Xw_mean, RA, Cc, ma);  %solve exterior orientation
%         err(3) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);
% 
%     end

    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor_l_pt(X3, Alph, Alph_T, norm_cw_tot, SUac_tot);
    [R,T] = epnp_lines.getrotT(XT_mean, RA_tot, Cc, m_tot);
    err(3) = epnp_lines.reprojection_error_usingRTAll(Xw, U, Xs,Xe,xs,xe,R,T,A);

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
if min(err)>THRESHOLD_REPROJECTION_ERROR %just compute if we do not have good solution in the previus cases
%     fprintf('solve 4\n');
    Km1=Km(:,end-3);
    Km2=Km(:,end-2);
    Km3=Km(:,end-1);
    Km4=Km(:,end);


    D=epnp_orig.compute_constraint_distance_orthog_4param_9eq_10unk(Km1,Km2,Km3,Km4);
    dsq=epnp_orig.define_distances_btw_control_points();
    lastcolumn=[-dsq',0,0,0]';
    D_=[D,lastcolumn];
    Kd=null(D_);

    P=epnp_orig.compute_permutation_constraint4(Kd);
    lambdas_=epnp_orig.kernel_noise(P,1);
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

%     if (length(Xw) <= 3)        
%         [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor(X4,Alph_S, SUac, norm_cw);
%         [R,T] = epnp_lines.getrotT(Xs_mean, RA, Cc, ms); %solve exterior orientation
%         err(4) = epnp_lines.reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A);
%     else
%         [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor(X4,Alph, SUac, norm_cw);
%         [R,T] = epnp_lines.getrotT(Xw_mean, RA, Cc, ma);  %solve exterior orientation
%         err(4) = epnp_lines.reprojection_error_usingRT(Xw,U,R,T,A);
%     end
    [Cc,Xc,sc] = epnp_lines.compute_norm_sign_scaling_factor_l_pt(X4, Alph, Alph_T, norm_cw_tot, SUac_tot);
    [R,T] = epnp_lines.getrotT(XT_mean, RA_tot, Cc, m_tot);
    err(4) = epnp_lines.reprojection_error_usingRTAll(Xw, U, Xs,Xe,xs,xe,R,T,A);


    sol(4).Xc=Xc;
    sol(4).Cc=Cc;
    sol(4).R=R;
    sol(4).T=T;
    sol(4).error=err(4);
    sol(4).betas=[beta1,beta2,beta3,beta4];
    sol(4).sc=sc;
    sol(4).Kernel=[Km1,Km2,Km3,Km4];
end
t2 = toc(t0);
t0 = tic;
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
% if (length(Xw) < 3)
%     [Xc_opt,R_opt,T_opt,err_opt,iter] = epnp_lines.optimize_betas_gauss_newton(Kernel,Cw,Beta0,Xs,Xe,xs,xe,A,Xw,U,Alph,Alph_S, Xs_mean, RA, ms);
% else
%     [Xc_opt,R_opt,T_opt,err_opt,iter] = epnp_lines.optimize_betas_gauss_newton(Kernel,Cw,Beta0,Xs,Xe,xs,xe,A,Xw,U,Alph, Alph_S, Xw_mean, RA, ma);
% end

[Xc_opt,R_opt,T_opt,err_opt,iter] = epnp_lines.optimize_betas_gauss_newton2(Kernel,Cw,Beta0,Xs,Xe,xs,xe,A,Xw,U, XT_mean, RA_tot, m_tot,Alph_T);

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

t3 = toc(t0);

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
