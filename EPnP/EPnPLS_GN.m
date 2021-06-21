function [R,t,ff, s] = EPnPLS_GN(Xw,xx,xs, xe, Xs, Xe, lnc)
    ff = 0;
    t0 = tic;
    
    if (nargin < 7)
        lnc = 1;
    end
    
    if (length(Xw(1,:)) + length(Xs(1,:)) < 6)
        R = [];
        t = [];
        return;
    end 

%preprocess with common shift
    meanT = mean([Xs Xe Xw], 2);
    Xs = Xs - repmat(meanT, 1, size(Xs, 2));
    Xe = Xe - repmat(meanT, 1, size(Xs, 2));
    Xw = Xw - repmat(meanT, 1, size(Xw, 2));
    
%normalize lines
    Xd = Xe - Xs;
    Xdn = sqrt(Xd(1, :).^2 + Xd(2, :).^2 + Xd(3, :).^2);
    Xd = Xd ./ repmat(Xdn, 3, 1);
    Xe = Xs + Xd;
        
    l2d = epnp_lines.pnl_preprocess(xs, xe);
    t1 = toc(t0);
    t0 = tic;
    [R,t,~,~,~,elt1,elt2,elt3,elt4]= epnp_lines.efficient_pnp_gauss_lines1(Xw, xx, Xs, Xe, xs, xe, diag([1 1 1]), l2d, lnc);  
    t2 = toc(t0);
    t0 = tic;
    
    Xsc = R*Xs + repmat(t, 1, size(Xs, 2));    
    Xec = R*Xe + repmat(t, 1, size(Xe, 2));      
    %1. find a point M : a middle point of a 3D line segment with endpoints
    %reprojected from image 2D endpoints
    %2. shift so that middle point of 3D line segment is iexactly same
    %point as M        

    [Xs, Xe] = epnp_lines.moveLineCloser2DVector(Xsc, Xec, xs, xe);

    c = -R'*t;
    Xs = R'*Xs + repmat(c, 1, size(Xs, 2));
    Xe = R'*Xe + repmat(c, 1, size(Xe, 2));
    t3 = toc(t0);
    t0 = tic;        
    [R,t]= epnp_lines.efficient_pnp_gauss_lines1(Xw, xx, Xs, Xe, xs, xe, diag([1 1 1]), l2d, lnc);        
    t4 = toc(t0);
        
    s = [t1,t2,t3,t4];
    %correct preprocess
    t = t - R*meanT;
    return    
end

function [Xs_new, Xe_new] = moveLineCloser2DNoR(Xsi, Xei, xs, xe, R, t)
    lineCoefs = cross(R*Xsi+t, R*Xei+t);
    lineCoefs = lineCoefs / norm(lineCoefs(1:2));    
    lineSpt = lineCoefs(1:2)*(-lineCoefs(3));
    lineDir = [-lineCoefs(2); lineCoefs(1)];
    ls = lineDir'*(xe-xs);
    if (ls<0)
        lineDir = -lineDir;
    end
    segLen = norm(xs-xe);
%     b = [xs-lineSpt; xe-lineSpt-lineDir*segLen];
%     A = [lineDir; lineDir];
%     sol = A\b;
    sol = lineDir'*(xs+xe-2*lineSpt)/2-segLen/2;
    xsn = lineSpt + sol*lineDir;
    xen = lineSpt + (sol+segLen)*lineDir;

    Xs_new = reproject2DtoLineNoR(xsn, Xsi, Xei, R, t);
    Xe_new = reproject2DtoLineNoR(xen, Xsi, Xei, R, t);
end


function [Xs_new, Xe_new, xsn, xen] = moveLineCloser2D(Xsi, Xei, xs, xe)
    lineCoefs = epnp_lines.cross1(Xsi, Xei);
    lineCoefs = lineCoefs / norm(lineCoefs(1:2));    
    lineSpt = lineCoefs(1:2)*(-lineCoefs(3));
    lineDir = [-lineCoefs(2); lineCoefs(1)];
    ls = lineDir'*(xe-xs);
    if (ls<0)
        lineDir = -lineDir;
    end
    segLen = norm(xs-xe);
%     b = [xs-lineSpt; xe-lineSpt-lineDir*segLen];
%     A = [lineDir; lineDir];
%     sol = A\b;
    sol = lineDir'*(xs+xe)/2-segLen/2;
    xsn = lineSpt + sol*lineDir;
    xen = lineSpt + (sol+segLen)*lineDir;
    
%     xm = 0.5*(xsn+xen);
%     Xm_new = reproject2DtoLine(xm, Xsi, Xei);
%     Xm = 0.5*(Xsi+Xei);
%     Shift = Xm_new - Xm;
%     Xs_new = Xsi+Shift;
%     Xe_new = Xei+Shift;
    Xs_new = reproject2DtoLine(xsn, Xsi, Xei);
    Xe_new = reproject2DtoLine(xen, Xsi, Xei);
end

function Xs_new = reproject2DtoLineNoR(xsi, Xsi, Xei, R, t)
    if (length(xsi) ~= size(xsi, 1))
        fprintf('ERROR in reproject2DtoLine \n');
    end
    if (length(xsi) == 2)
        xsi = [xsi; 1];
    end
    Xdir = R*(Xei-Xsi);
    A = [xsi -Xdir];
    b = R*Xsi + t;
    sol = A\b;
    Xs_new = R*Xsi + t + Xdir*sol(2);
end

function Xs_new = reproject2DtoLine(xsi, Xsi, Xei)
    if (length(xsi) ~= size(xsi, 1))
        fprintf('ERROR in reproject2DtoLine \n');
    end
    if (length(xsi) == 2)
        xsi = [xsi; 1];
    end
    Xdir = Xei-Xsi;
    A = [xsi -Xdir];
    b = Xsi;
    sol = A\b;
    Xs_new = Xsi + Xdir*sol(2);
end

function Xs_new = reproject2DtoLineIn2D(xsi, Xsi, Xei)
    if (length(xsi) ~= size(xsi, 1))
        fprintf('ERROR in reproject2DtoLine \n');
    end
    lineCoefs = cross(Xsi, Xei);
    lineCoefs = lineCoefs / norm(lineCoefs(1:2));    
    lineSpt = lineCoefs(1:2)*(-lineCoefs(3));
    lineDir = [-lineCoefs(2); lineCoefs(1)];
    A = -lineDir;
    b = lineSpt - xsi;
    sol = A\b;
    linePt = lineSpt + lineDir*sol(1);
    Xs_new = reproject2DtoLine(linePt, Xsi, Xei);
end

function [Xs, Xe] = shiftLinesToPoint(pt, Xs, Xe)
    %1) determine closest point on each line to the point pt
    %2) shift the line so that middle point of segment is there
    for i = 1:size(Xs, 2)
        Xdir = Xe(:, i) - Xs(:, i);
        sol = Xdir\(pt - Xs(:, i));
        Cpt = Xs(:, i)+Xdir*sol;
        Mpt = 0.5*(Xe(:, i)+Xs(:, i));
        Shift = Mpt - Cpt;
        Xe(:, i) = Xe(:, i)-Shift;
        Xs(:, i) = Xs(:, i)-Shift;
    end
end

