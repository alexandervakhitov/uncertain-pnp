function [R, T, inliers, trialcount] = ransac(K, X, x, fittingfn, s, t, ...
    scales_2d, idsampl, maxTrials)

    % Test number of parameters
    narginchk ( 6, 9 ) ;
    
    [~, npts] = size(x);
     
    if nargin < 8; idsampl = 1:npts;    end;
    if nargin < 9; maxTrials = 1000;    end;
    
    p = 0.99;         % Desired probability of choosing at least one sample
                      % free from outliers

    bestR = NaN;      % Sentinel value allowing detection of solution failure.
    bestT = NaN;
    trialcount = 0;
    bestscore =  0;
    
    %we assume 50% of outliers
    %N = 1;            % Dummy initialisation for number of trials.
    
    N = log(1-p)/log(1-0.5^s);
    nsampl = length(idsampl);
    
    while N > trialcount
      
        % Generate s random indicies in the range 1..npts
        ind = randsample(nsampl, s);
        ind = idsampl(ind);
        
        try
            %K = [f 0 0;0 f 0; 0 0 1];
            u = K \ [x(:,ind);ones(1,numel(ind))];
            u = u(1:2,:);
            [R, T] = fittingfn(X(:,ind), u);
        catch
            R = [];
            T = [];
        end
        [inliers, R, T] = ransac_eval(K, R, T, X, x, t, scales_2d);

        ninliers = length(inliers);
        
        if ninliers > bestscore    % Largest set of inliers so far...
            bestscore = ninliers;  % Record data for this model
            bestinliers = inliers;
            bestR = R;
            bestT = T;
            
            % Update estimate of N, the number of trials to ensure we pick,
            % with probability p, a data set with no outliers.
            fracinliers =  ninliers/nsampl;
            pNoOutliers = 1 -  fracinliers^s;
            pNoOutliers = max(eps, pNoOutliers);  % Avoid division by -Inf
            pNoOutliers = min(1-eps, pNoOutliers);% Avoid division by 0.
            N = log(1-p)/log(pNoOutliers);
        end
        
        trialcount = trialcount+1;
        
        
        % Safeguard against being stuck in this loop forever
        if trialcount > maxTrials
            break
        end
    end
    
    if ~isnan(bestR)   % We got a solution
        R = bestR;
        T = bestT;
        inliers = bestinliers;
    else
        R = [];
        T = [];
        inliers = [];
    end
end

