function [err, U_S_rep, U_E_rep, errs] = reprojection_error_usingRTAll(Xw, U, Xs,Xe,l2d,R,T,A)
   U_S_rep = epnp_lines.projectPts(Xs, A, R, T);
   U_E_rep = epnp_lines.projectPts(Xe, A, R, T); 
   
   err_s = l2d(1:2, :)' * U_S_rep + l2d(3, :)';
   err_e = l2d(1:2, :)' * U_E_rep + l2d(3, :)';
   err_l = err_s.^2 + err_e.^2;   
   
   Urep = epnp_lines.projectPts(Xw, eye(3), R, T);    
   err_ = (U(1,:)-Urep(1,:)).^2+(U(2,:)-Urep(2,:)).^2;

   err = sqrt((sum(err_l) + sum(err_))/(size(Xs, 2)+size(Xw, 2)));
end