function [err, U_S_rep, U_E_rep, errs] = reprojection_error_usingRTLines(Xs,Xe,xs,xe,R,T,A)
   U_S_rep = epnp_lines.projectPts(Xs, A, R, T);
   U_E_rep = epnp_lines.projectPts(Xe, A, R, T); 
   Xsc = A*R*Xs + repmat(A*T, 1, size(xs, 2));
   Xec = A*R*Xe + repmat(A*T, 1, size(xe, 2));
   l2dP = cross(Xsc, Xec);
   l2dn = sqrt(l2dP(1,:).^2 + l2dP(2,:).^2);
   l2dP1 = l2dP./repmat(l2dn, 3, 1);
   prerrs_xs = l2dP1(1,:).*xs(1,:) + l2dP1(2,:).*xs(2,:) + l2dP1(3, :);
   prerrs_xe = l2dP1(1,:).*xe(1,:) + l2dP1(2,:).*xe(2,:) + l2dP1(3, :);
   errs = prerrs_xs.^2+prerrs_xe.^2;   
   err = sum(errs)/size(Xs, 2);
end