function Urep = projectPts(Xw, A, R, T)
    %clear all; close all; load reprojection_error_usingRT;
    n=size(Xw,2);

    PR=A*[R];
    PT = A*T;
%     Xw_h=[Xw; ones(1,n)];
    Urep_ = PR*Xw + repmat(PT, 1, n);

    %project reference points into the image plane
    Urep = Urep_(1:2, :)./repmat(Urep_(3, :), 2, 1);
%     Urep(1,:)=Urep_(1,:)./Urep_(3,:);
%     Urep(2,:)=Urep_(2,:)./Urep_(3,:);

end