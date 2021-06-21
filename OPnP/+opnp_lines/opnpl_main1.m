function [R0, t0, fail_flag] = opnpl_main1(U,u,Xs, Xe, l2d, label_polish)
    %will polish the solution in default
    if nargin < 6
        label_polish = 'polish';
    end    
    n = size(U,2);
    nl = size(Xs, 2);
    %homogeneous coordinate
    if size(u,1) > 2
        u = u(1:2,:);
    end

    %3D points after translation to centroid
    Ucent = mean([U Xs Xe],2);
    Um = U - repmat(Ucent,1,n);

    xm = Um(1,:)'; ym = Um(2,:)'; zm = Um(3,:)';
    x = U(1,:)'; y = U(2,:)'; z = U(3,:)';
    u1 = u(1,:)'; v1 = u(2,:)';


    %construct matrix N: 2n*11
    M = zeros(2*n,11);
    %1 a2 ab ac ad b2 bc bd c2 cd d2
    M(1:2:end,:) = [ u1, u1.*zm - x, 2*u1.*ym, - 2*z - 2*u1.*xm, 2*y, - x - u1.*zm, -2*y, 2*u1.*xm - 2*z, x - u1.*zm, 2*u1.*ym, x + u1.*zm];
    M(2:2:end,:) = [ v1, v1.*zm - y, 2*z + 2*v1.*ym, -2*v1.*xm, -2*x, y - v1.*zm, -2*x, 2*v1.*xm, - y - v1.*zm, 2*v1.*ym - 2*z, y + v1.*zm];
    N = zeros(2*n, 2);
    N(1:2:end, 1:2) = [-ones(n, 1), zeros(n, 1)];
    N(2:2:end, 1:2) = [zeros(n, 1), -ones(n, 1)];
    
%     MTN = [sum(N(1:2:end,:)); sum(N(2:2:end,:))];
    

    %construct matrix Q: 11*11
%     Q = N'*N - 1/n*(MTN')*MTN;
    %
    RtoPol = opnp_lines.generateRToPol();

    L = l2d';        

    Ml = zeros(2*nl, 11);
    %1 a2 ab ac ad b2 bc bd c2 cd d2
    x = Xs(1, :)';
    y = Xs(2, :)';
    z = Xs(3, :)';
    xm = Ucent(1);
    ym = Ucent(2);
    zm = Ucent(3);
    l1 = l2d(1,:)';
    l2 = l2d(2,:)';
    l3 = l2d(3,:)';
    Ml(1:nl,:) = [ l3, l1.*x + l2.*y + l3.*z - l3.*zm, -2*l2.*z+2*l3.*y-2*l3.*ym, 2*l1.*z - 2*l3.*x + 2*l3.*xm, -2*l1.*y+2*l2.*x, l1.*x-l2.*y - l3.*z + l3.*zm, 2*l1.*y+2*l2.*x, 2*l1.*z + 2*l3.*x - 2*l3.*xm, -l1.*x+l2.*y - l3.*z + l3.*zm, 2*l2.*z+2*l3.*y-2*l3.*ym, -l1.*x-l2.*y+l3.*z - l3.*zm];
    x = Xe(1, :)';
    y = Xe(2, :)';
    z = Xe(3, :)';
    xm = Ucent(1);
    ym = Ucent(2);
    zm = Ucent(3);    
    Ml(nl+1:2*nl,:) = [ l3, l1.*x + l2.*y + l3.*z - l3.*zm, -2*l2.*z+2*l3.*y-2*l3.*ym, 2*l1.*z - 2*l3.*x + 2*l3.*xm, -2*l1.*y+2*l2.*x, l1.*x-l2.*y - l3.*z + l3.*zm, 2*l1.*y+2*l2.*x, 2*l1.*z + 2*l3.*x - 2*l3.*xm, -l1.*x+l2.*y - l3.*z + l3.*zm, 2*l2.*z+2*l3.*y-2*l3.*ym, -l1.*x-l2.*y+l3.*z - l3.*zm];

    
%     Lrep = [repmat(L(:, 1), 1, 3) repmat(L(:, 2), 1, 3) repmat(L(:, 3), 1, 3)];
%     N011 = repmat(Xs', 1, 3).*Lrep;
%     N012 = repmat(Xe', 1, 3).*Lrep;    
    Lbig = repmat(L(:,1:2), 2, 1);
%     Ml = N0*RtoPol';
    Nl = Lbig;

    Ntot = [N; Nl];
    Mtot = [M; Ml];

    Ntt = Ntot'*Ntot;
    if (sum(isnan(Ntt(:))) || rank(Ntt) < 2)
        R0 = [];
        t0 = [];
        fail_flag = 1;
        return;
    end
    Ptot = -pinv(Ntt)*(Ntot'*Mtot);

    Q = Mtot'*Mtot + Mtot'*Ntot*Ptot;

    const = Q(1,1);
    q = Q(1,2:end);
    Q = Q(2:end,2:end);

    Q11 = Q(1,1); Q12 = Q(1,2); Q13 = Q(1,3); Q14 = Q(1,4); Q15 = Q(1,5); Q16 = Q(1,6); Q17 = Q(1,7); Q18 = Q(1,8); Q19 = Q(1,9); Q110 = Q(1,10); 
                  Q22 = Q(2,2); Q23 = Q(2,3); Q24 = Q(2,4); Q25 = Q(2,5); Q26 = Q(2,6); Q27 = Q(2,7); Q28 = Q(2,8); Q29 = Q(2,9); Q210 = Q(2,10); 
                                Q33 = Q(3,3); Q34 = Q(3,4); Q35 = Q(3,5); Q36 = Q(3,6); Q37 = Q(3,7); Q38 = Q(3,8); Q39 = Q(3,9); Q310 = Q(3,10); 
                                              Q44 = Q(4,4); Q45 = Q(4,5); Q46 = Q(4,6); Q47 = Q(4,7); Q48 = Q(4,8); Q49 = Q(4,9); Q410 = Q(4,10); 
                                                            Q55 = Q(5,5); Q56 = Q(5,6); Q57 = Q(5,7); Q58 = Q(5,8); Q59 = Q(5,9); Q510 = Q(5,10);
                                                                          Q66 = Q(6,6); Q67 = Q(6,7); Q68 = Q(6,8); Q69 = Q(6,9); Q610 = Q(6,10);
                                                                                        Q77 = Q(7,7); Q78 = Q(7,8); Q79 = Q(7,9); Q710 = Q(7,10);
                                                                                                      Q88 = Q(8,8); Q89 = Q(8,9); Q810 = Q(8,10);
                                                                                                                    Q99 = Q(9,9); Q910 = Q(9,10);
                                                                                                                                  Q1010 = Q(10,10);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6); q7 = q(7); q8 = q(8); q9 = q(9); q10 = q(10);

    %variable sequence
    %[ a^3, a^2*b, a^2*c, a^2*d, a*b^2, a*b*c, a*b*d, a*c^2, a*c*d, a*d^2, a, b^3, b^2*c, b^2*d, b*c^2, b*c*d, b*d^2, b, c^3, c^2*d, c*d^2, c, d^3, d]

    coeff1 = [ 4*Q11, 6*Q12, 6*Q13, 6*Q14, 4*Q15 + 2*Q22, 4*Q16 + 4*Q23, 4*Q17 + 4*Q24, 4*Q18 + 2*Q33, 4*Q19 + 4*Q34, 4*Q110 + 2*Q44, 4*q1, 2*Q25, 2*Q26 + 2*Q35, 2*Q27 + 2*Q45, 2*Q28 + 2*Q36, 2*Q29 + 2*Q37 + 2*Q46, 2*Q210 + 2*Q47, 2*q2, 2*Q38, 2*Q39 + 2*Q48, 2*Q310 + 2*Q49, 2*q3, 2*Q410, 2*q4];
    coeff2 = [ 2*Q12, 4*Q15 + 2*Q22, 2*Q16 + 2*Q23, 2*Q17 + 2*Q24, 6*Q25, 4*Q26 + 4*Q35, 4*Q27 + 4*Q45, 2*Q28 + 2*Q36, 2*Q29 + 2*Q37 + 2*Q46, 2*Q210 + 2*Q47, 2*q2, 4*Q55, 6*Q56, 6*Q57, 4*Q58 + 2*Q66, 4*Q59 + 4*Q67, 4*Q510 + 2*Q77, 4*q5, 2*Q68, 2*Q69 + 2*Q78, 2*Q610 + 2*Q79, 2*q6, 2*Q710, 2*q7];
    coeff3 = [ 2*Q13, 2*Q16 + 2*Q23, 4*Q18 + 2*Q33, 2*Q19 + 2*Q34, 2*Q26 + 2*Q35, 4*Q28 + 4*Q36, 2*Q29 + 2*Q37 + 2*Q46, 6*Q38, 4*Q39 + 4*Q48, 2*Q310 + 2*Q49, 2*q3, 2*Q56, 4*Q58 + 2*Q66, 2*Q59 + 2*Q67, 6*Q68, 4*Q69 + 4*Q78, 2*Q610 + 2*Q79, 2*q6, 4*Q88, 6*Q89, 4*Q810 + 2*Q99, 4*q8, 2*Q910, 2*q9];
    coeff4 = [ 2*Q14, 2*Q17 + 2*Q24, 2*Q19 + 2*Q34, 4*Q110 + 2*Q44, 2*Q27 + 2*Q45, 2*Q29 + 2*Q37 + 2*Q46, 4*Q210 + 4*Q47, 2*Q39 + 2*Q48, 4*Q310 + 4*Q49, 6*Q410, 2*q4, 2*Q57, 2*Q59 + 2*Q67, 4*Q510 + 2*Q77, 2*Q69 + 2*Q78, 4*Q610 + 4*Q79, 6*Q710, 2*q7, 2*Q89, 4*Q810 + 2*Q99, 6*Q910, 2*q9, 4*Q1010, 4*q10];

    %normalize 
    coeff1 = coeff1/norm(coeff1); 
    coeff2 = coeff2/norm(coeff2); 
    coeff3 = coeff3/norm(coeff3); 
    coeff4 = coeff4/norm(coeff4); 

    %call grobner basis solver
    [xx yy zz tt, fail_flag] = opnp_orig.GB_Solver_3Order_4Variable_Symmetry(coeff1, coeff2, coeff3, coeff4);
    if (fail_flag == 1)
        R0 = [];
        t0 = [];
        fail_flag = 1;
        return;
    end

    %filter out repetitive solutions
    [xx yy zz tt] = opnp_orig.PnP_FilterOutRepetitive(xx,yy,zz,tt);

%     polish solutions
    if strcmp(label_polish,'polish')
        for jj = 1:length(xx)
            [xx(jj) yy(jj) zz(jj) tt(jj) label_psd(jj)] = opnp_orig.PnP_Polish(Q,q,[xx(jj) yy(jj) zz(jj) tt(jj)]);
        end
    end

    index = 1;
    error = [];
    %recover R and t
    for i = 1:length(xx)    
        a = xx(i); b = yy(i); c = zz(i); d = tt(i); 

        lambda1 = 1/(a^2+b^2+c^2+d^2);
        if lambda1 > 1e10
            continue;
        end

        vec = [1 a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2];        

        a = a*sqrt(lambda1); b = b*sqrt(lambda1); c = c*sqrt(lambda1); d = d*sqrt(lambda1);

        R = [a^2+b^2-c^2-d^2     2*b*c-2*a*d     2*b*d+2*a*c
             2*b*c+2*a*d         a^2-b^2+c^2-d^2 2*c*d-2*a*b
             2*b*d-2*a*c         2*c*d+2*a*b   a^2-b^2-c^2+d^2];

        t = [Ptot*vec'*lambda1; lambda1 - R(3, :)*Ucent];

        proj = R*U+t*ones(1,n);

        if min(proj(3,:)) < 0
            continue;
        end

        %proj = proj./repmat(proj(3,:),3,1);    
        %err_proj = norm(proj(1:2,:)-u(1:2,:),'fro')/n;
        err = vec*[const q;q.' Q]*vec.';

        Rot(:,:,index) = R; 
        trans(:,index) = t; 
        error(index) = err; %error1(index) = err_proj;
        index = index + 1;
    end

    if isempty(error)
        R0 = []; t0 = []; error0 = []; flag = 1;
        return;
    end

    error0 = min(abs(error));

    %the minimial-case or the slightly over-constrained cases
    if  n <= 6
        R0 = Rot;
        t0 = trans;
        flag = 0;
    %in the noise and overconstrained case
    else
        rr = find(error <= error0); 
        R0 = Rot(:,:,rr);
        t0 = trans(:,rr);
        flag = 0;
    end


end