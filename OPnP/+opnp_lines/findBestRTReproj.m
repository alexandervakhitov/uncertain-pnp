function [minInd Xsf Xef Lcf] = findBestRTReproj(R0, t0, xs, xe, Xs, Xe, U, u)
    if (size(R0, 1) == 0)
        minInd = 0;
        Xsf = [];
        Xef = [];
        Lcf = [];
        return;        
    end
    minErr = 1e10;
    minInd = 1;
    for i = 1:size(R0, 3)
        R = R0(:, :, i);
        t = t0(:, i); 
        
        if (size(U, 2) > 0)
            projs = R*U + repmat(t, 1, size(U, 2));
            projs(1:2, :) = projs(1:2, :) ./ repmat(projs(3, :), 2, 1);
            res = u(1:2, :) - projs(1:2, :);
            aggPtErr = mean(sqrt(res(1, :).^2 + res(2, :).^2));
        else
            aggPtErr = 0;
        end
        
        Xsp = R*Xs + repmat(t, 1, size(Xs, 2));
        Xep = R*Xe + repmat(t, 1, size(Xe, 2));        
        
        LcUnNormed = cross(Xsp, Xep);
        LcNorms = sqrt(LcUnNormed(1,:).^2 + LcUnNormed(2,:).^2);
        LcNormed = LcUnNormed ./ repmat(LcNorms, 3, 1);
        
        if (size(R0, 3) == 1)
            minInd = 1;
            Xsf = Xsp;
            Xef = Xep;
            Lcf = LcNormed;
            return;
        end
        
        aggErrLi = sum(abs(LcNormed(1,:).*xs(1,:) + LcNormed(2,:).*xs(2,:) + LcNormed(3,:)));
        aggErrLi = aggErrLi + sum(abs(LcNormed(1,:).*xe(1,:) + LcNormed(2,:).*xe(2,:) + LcNormed(3,:)));
        errVal = aggPtErr + aggErrLi;
        if (errVal < minErr)
            minErr = errVal;
            minInd = i;
            Xsf = Xsp;
            Xef = Xep;
            Lcf = LcNormed;
        end
    end    
end
