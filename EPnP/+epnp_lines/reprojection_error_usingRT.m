function [err,Urep]=reprojection_error_usingRT(Xw,U,R,T,A)

    Urep = epnp_lines.projectPts(Xw, A, R, T);

    %reprojection error
    err_=sqrt((U(1,:)-Urep(1,:)).^2+(U(2,:)-Urep(2,:)).^2);
    n = length(err_);
    err=sum(err_)/n;

end
