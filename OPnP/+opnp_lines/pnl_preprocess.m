function l2d = pnl_preprocess(xs, xe)
    if (size(xs, 1) == 2)
        xs = [xs; ones(1, size(xs, 2))];        
    end
    if (size(xe, 1) == 2)
        xe = [xe; ones(1, size(xe, 2))];        
    end    
    
    l2d = cross(xs, xe);
    l2dn = sqrt(l2d(1,:).^2+l2d(2,:).^2);
    l2d = l2d ./ repmat(l2dn, 3, 1);
        
end