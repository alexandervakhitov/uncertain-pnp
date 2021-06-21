function [Xsn, Xen] = moveLineCloser2DLcVector(Xsi, Xei, xs, xe, Lc)
    LineSpt = -Lc(1:2, :).*repmat(Lc(3, :), 2, 1);
    Ld = [-Lc(2, :); Lc(1, :)];
    xd = xe - xs;

    ls = sign(Ld(1,:).*xd(1,:) + Ld(2, :).*xd(2, :));
    Ld = Ld.*repmat(ls, 2, 1);

    SegLen = sqrt(xd(1,:).^2 + xd(2,:).^2);
    xsc = (xs+xe)/2;
    p = Ld(1,:).*xsc(1,:) + Ld(2, :).*xsc(2, :);
    Sol = p-SegLen/2;
    Sol2 = repmat(Sol, 2, 1);
    xsn = LineSpt + Ld.*Sol2;
    Sol2Seg = repmat(Sol+SegLen, 2, 1);
    xen = LineSpt + Ld.*Sol2Seg;

    Xd = Xei-Xsi;
    a1 = [Xsi(1, :); Xd(1, :); -xsn(1,:)];
    a2 = [Xsi(2, :); Xd(2, :); -xsn(2,:)];
    p = cross(a1, a2);
    q = p(2,:)./p(1,:);
    Xsn = repmat(q,3,1).*Xd + Xsi;

    a1 = [Xsi(1, :); Xd(1, :); -xen(1,:)];
    a2 = [Xsi(2, :); Xd(2, :); -xen(2,:)];
    p = cross(a1, a2);
    q = p(2,:)./p(1,:);
    Xen = repmat(q,3,1).*Xd + Xsi;


end
