function [R0 t0 error0 flag] = OPnPL(U,u,xs, xe, Xs, Xe,label_polish)
%will polish the solution in default
    if nargin < 7
        label_polish = 'polish';
    end    
    n = size(U,2);
    
    %normalize lines
    Xd = Xe - Xs;
    Xdn = sqrt(Xd(1, :).^2 + Xd(2, :).^2 + Xd(3, :).^2);
    Xd = Xd ./ repmat(Xdn, 3, 1);
    Xe = Xs + Xd;
    
    %check for planar case
    AllPts = [Xs Xe U];
    AllPts_cent = AllPts - repmat(mean(AllPts, 2), 1, size(AllPts, 2)); 
    PtCorr = AllPts_cent*AllPts_cent';
    [PU, PS, PV] = svd(PtCorr*PtCorr');        
    if (abs(PS(3,3)) < 1e-10)
        isPlanar = 1;
    else
        isPlanar = 0;
    end
    
    l2d = opnp_lines.pnl_preprocess(xs, xe);

    if (n>2)
        [R0 t0] = OPnP(U,u);
        error0 = 0;
        flag = 0;
        as = zeros(size(R0, 3), 1);
    else
        if (isPlanar)
            [R0 t0 error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, 2);
            while (size(R0, 3) == 0)
                [R0 t0 error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, 1);
            end            
        else
            [R0 t0 error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d);
            while (size(R0, 3) == 0)
                [R0 t0 error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, 1);
            end
        end       
    end   
    [minInd Xsc Xec Lc] = opnp_lines.findBestRTReproj(R0, t0, xs, xe, Xs, Xe, U, u);
    [Xs1, Xe1] = opnp_lines.moveLineCloser2DLcVector(Xsc, Xec, xs, xe, Lc);
    R = R0(:, :, minInd);
    t = t0(:, minInd);    
    Xs = R'*Xs1 - repmat(R'*t, 1, size(Xs1, 2));
    Xe = R'*Xe1 - repmat(R'*t, 1, size(Xs1, 2));  

    if (n > 0)
        [R0 t0] = opnp_lines.opnpl_main(U, u, Xs, Xe, l2d, label_polish);            
        error0 = 0;
        flag = 0;            
    else               
        if (isPlanar)
            [R0n t0n error0n flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, 2);
            if (error0n > error0)
                R0n = R0;
                t0n = t0;
            end
        else
            [R0n t0n error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, as(minInd));
            while (size(R0n, 3) == 0)
                [R0n t0n error0 flag as] = opnp_lines.opnl_main(Xs, Xe, l2d, 1);
            end
        end   
        if (size(R0n, 3) > 0)
            R0 = R0n;
            t0 = t0n;
        end
    end      
    
%     check for planar case    
    if (isPlanar);
        PS;
        R0o = R0;
        t0o = t0;
        initSolNum = size(R0, 3);
        R0 = zeros(3, 3, 2*initSolNum);
        t0 = zeros(3, 2*initSolNum);
        for i = 1:initSolNum
            Ri = R0o(:, :, i);
            ti = t0o(:, i);
            Ri = Ri*PV;            
            R0(:, :, 2*i-1) = Ri*PV';
            t0(:, 2*i-1) = ti;
            R0(:, :, 2*i) = [-Ri(:, 1) -Ri(:, 2) Ri(:, 3)]*PV';
            t0(:, 2*i) = -ti;
        end
    end
end
function optInd = findBestRT(R01, R, t01, t)
    minDiff = 1e10;
    optInd = 1;
    for jj = 1:size(R01, 3)
        diff = norm(R01(:,:, jj) - R) + norm(t01(:, jj) - t);
        if (diff < minDiff)
            minDiff = diff;
            optInd = jj;
        end
    end
end

function [Xs_new, Xe_new] = moveLineCloser2D(Xsi, Xei, xs, xe)
    lineCoefs = cross(Xsi, Xei);
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
    
    Xs_new = reproject2DtoLine(xsn, Xsi, Xei);
    Xe_new = reproject2DtoLine(xen, Xsi, Xei);
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

    
      
