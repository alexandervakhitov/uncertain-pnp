function [xx yy zz tt label_psd obj_cur] = PnP_Polish(Q,w,q0)
%q=[a b c d];
%vec = [1 a^2 ab ac ad b^2 bc bd c^2 cd d^2];
%objective function: vec.'*Q*vec
%refine the solution by using damped Newton method

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
q1 = w(1); q2 = w(2); q3 = w(3); q4 = w(4); q5 = w(5); q6 = w(6); q7 = w(7); q8 = w(8); q9 = w(9); q10 = w(10);

%maximum allowed iterations (one-step in our implementation)
maxItr = 1; 
%damped factor
lambda = 1e-8;
%max lambda
maxLambda = 1e5;

%min lambda
minLambda = 1e-8;

%flag 
flag = 0; 

q = q0;
a = q0(1); b = q0(2); c = q0(3); d = q0(4);
vec = [a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2].';
obj_pre = vec.'*Q*vec + 2*w*vec;

itr = 1; 
%iteration 
while itr <= maxItr 
        
    a = q(1); b = q(2); c = q(3); d = q(4);
    vec = [a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2 1].';
    vec3 = [a^3, a^2*b, a^2*c, a^2*d, a*b^2, a*b*c, a*b*d, a*c^2, a*c*d, a*d^2, a, b^3, b^2*c, b^2*d, b*c^2, b*c*d, b*d^2, b, c^3, c^2*d, c*d^2, c, d^3, d].';

    %coefficients for gradient
    coeff1 = [ 4*Q11, 6*Q12, 6*Q13, 6*Q14, 4*Q15 + 2*Q22, 4*Q16 + 4*Q23, 4*Q17 + 4*Q24, 4*Q18 + 2*Q33, 4*Q19 + 4*Q34, 4*Q110 + 2*Q44, 4*q1, 2*Q25, 2*Q26 + 2*Q35, 2*Q27 + 2*Q45, 2*Q28 + 2*Q36, 2*Q29 + 2*Q37 + 2*Q46, 2*Q210 + 2*Q47, 2*q2, 2*Q38, 2*Q39 + 2*Q48, 2*Q310 + 2*Q49, 2*q3, 2*Q410, 2*q4];
    coeff2 = [ 2*Q12, 4*Q15 + 2*Q22, 2*Q16 + 2*Q23, 2*Q17 + 2*Q24, 6*Q25, 4*Q26 + 4*Q35, 4*Q27 + 4*Q45, 2*Q28 + 2*Q36, 2*Q29 + 2*Q37 + 2*Q46, 2*Q210 + 2*Q47, 2*q2, 4*Q55, 6*Q56, 6*Q57, 4*Q58 + 2*Q66, 4*Q59 + 4*Q67, 4*Q510 + 2*Q77, 4*q5, 2*Q68, 2*Q69 + 2*Q78, 2*Q610 + 2*Q79, 2*q6, 2*Q710, 2*q7];
    coeff3 = [ 2*Q13, 2*Q16 + 2*Q23, 4*Q18 + 2*Q33, 2*Q19 + 2*Q34, 2*Q26 + 2*Q35, 4*Q28 + 4*Q36, 2*Q29 + 2*Q37 + 2*Q46, 6*Q38, 4*Q39 + 4*Q48, 2*Q310 + 2*Q49, 2*q3, 2*Q56, 4*Q58 + 2*Q66, 2*Q59 + 2*Q67, 6*Q68, 4*Q69 + 4*Q78, 2*Q610 + 2*Q79, 2*q6, 4*Q88, 6*Q89, 4*Q810 + 2*Q99, 4*q8, 2*Q910, 2*q9];
    coeff4 = [ 2*Q14, 2*Q17 + 2*Q24, 2*Q19 + 2*Q34, 4*Q110 + 2*Q44, 2*Q27 + 2*Q45, 2*Q29 + 2*Q37 + 2*Q46, 4*Q210 + 4*Q47, 2*Q39 + 2*Q48, 4*Q310 + 4*Q49, 6*Q410, 2*q4, 2*Q57, 2*Q59 + 2*Q67, 4*Q510 + 2*Q77, 2*Q69 + 2*Q78, 4*Q610 + 4*Q79, 6*Q710, 2*q7, 2*Q89, 4*Q810 + 2*Q99, 6*Q910, 2*q9, 4*Q1010, 4*q10];
    
    %gradient
    g = [coeff1*vec3; coeff2*vec3; coeff3*vec3; coeff4*vec3];
  
    %coefficients for hessian
    h11 = [ 12*Q11, 12*Q12, 12*Q13, 12*Q14, 4*Q15 + 2*Q22, 4*Q16 + 4*Q23, 4*Q17 + 4*Q24, 4*Q18 + 2*Q33, 4*Q19 + 4*Q34, 4*Q110 + 2*Q44, 4*q1]; 
    h12 = [ 6*Q12, 8*Q15 + 4*Q22, 4*Q16 + 4*Q23, 4*Q17 + 4*Q24, 6*Q25, 4*Q26 + 4*Q35, 4*Q27 + 4*Q45, 2*Q28 + 2*Q36, 2*Q29 + 2*Q37 + 2*Q46, 2*Q210 + 2*Q47, 2*q2]; 
    h13 = [ 6*Q13, 4*Q16 + 4*Q23, 8*Q18 + 4*Q33, 4*Q19 + 4*Q34, 2*Q26 + 2*Q35, 4*Q28 + 4*Q36, 2*Q29 + 2*Q37 + 2*Q46, 6*Q38, 4*Q39 + 4*Q48, 2*Q310 + 2*Q49, 2*q3]; 
    h14 = [ 6*Q14, 4*Q17 + 4*Q24, 4*Q19 + 4*Q34, 8*Q110 + 4*Q44, 2*Q27 + 2*Q45, 2*Q29 + 2*Q37 + 2*Q46, 4*Q210 + 4*Q47, 2*Q39 + 2*Q48, 4*Q310 + 4*Q49, 6*Q410, 2*q4];
 
    h22 = [ 4*Q15 + 2*Q22, 12*Q25, 4*Q26 + 4*Q35, 4*Q27 + 4*Q45, 12*Q55, 12*Q56, 12*Q57, 4*Q58 + 2*Q66, 4*Q59 + 4*Q67, 4*Q510 + 2*Q77, 4*q5];
    h23 = [ 2*Q16 + 2*Q23, 4*Q26 + 4*Q35, 4*Q28 + 4*Q36, 2*Q29 + 2*Q37 + 2*Q46, 6*Q56, 8*Q58 + 4*Q66, 4*Q59 + 4*Q67, 6*Q68, 4*Q69 + 4*Q78, 2*Q610 + 2*Q79, 2*q6];
    h24 = [ 2*Q17 + 2*Q24, 4*Q27 + 4*Q45, 2*Q29 + 2*Q37 + 2*Q46, 4*Q210 + 4*Q47, 6*Q57, 4*Q59 + 4*Q67, 8*Q510 + 4*Q77, 2*Q69 + 2*Q78, 4*Q610 + 4*Q79, 6*Q710, 2*q7];
    
    h33 = [ 4*Q18 + 2*Q33, 4*Q28 + 4*Q36, 12*Q38, 4*Q39 + 4*Q48, 4*Q58 + 2*Q66, 12*Q68, 4*Q69 + 4*Q78, 12*Q88, 12*Q89, 4*Q810 + 2*Q99, 4*q8]; 
    h34 = [ 2*Q19 + 2*Q34, 2*Q29 + 2*Q37 + 2*Q46, 4*Q39 + 4*Q48, 4*Q310 + 4*Q49, 2*Q59 + 2*Q67, 4*Q69 + 4*Q78, 4*Q610 + 4*Q79, 6*Q89, 8*Q810 + 4*Q99, 6*Q910, 2*q9];
 
    h44 = [ 4*Q110 + 2*Q44, 4*Q210 + 4*Q47, 4*Q310 + 4*Q49, 12*Q410, 4*Q510 + 2*Q77, 4*Q610 + 4*Q79, 12*Q710, 4*Q810 + 2*Q99, 12*Q910, 12*Q1010, 4*q10];
 
    %hessian matrix
    H = [h11*vec h12*vec h13*vec h14*vec
         h12*vec h22*vec h23*vec h24*vec
         h13*vec h23*vec h33*vec h34*vec
         h14*vec h24*vec h34*vec h44*vec];
    
    %not positive definite Hessian
    eH = abs(eig(H));
    if min(eH)/max(eH) < -1e-3
        label_psd = 0; 
        xx = 0; yy = 0; zz = 0; tt = 0;
        obj_cur = inf;
        return;
    else
        label_psd = 1;
    end
    
    q_temp = q;
    while (lambda < maxLambda)
        %increment
        delta = (H+lambda*eye(4))\g;

        %update parameter
        q = q_temp - delta';
         
        a = q(1); b = q(2); c = q(3); d = q(4);
        vec = [a^2 a*b a*c a*d b^2 b*c b*d c^2 c*d d^2].';
        
        %evaluate the sampson error 
        obj_cur = vec.'*Q*vec + 2*w*vec;

        %check convergence
        if obj_cur >= obj_pre
            lambda = 10*lambda;
            continue;
        else
            obj_pre = obj_cur; 
            lambda = 0.1*lambda;
            break;
        end        
    end
    if lambda >= maxLambda
        q = q_temp;
        break;
    end
    
    if lambda <= minLambda
        lambda = minLambda;
    end
    itr = itr + 1;
end

xx = q(1); yy = q(2); zz = q(3); tt = q(4);
end

