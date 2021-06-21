function l2d = computeImageLines(xs, xe)
    %only for u0 = v0 = 0
    nl = size(xs, 2);
    l2d = zeros(3, nl);
    for i = 1:nl
        lineCoefs = epnp_lines.cross1([xs(:, i); 1], [xe(:, i); 1]);
        lineCoefs = lineCoefs / norm(lineCoefs(1:2));
        l2d(:, i) = lineCoefs(:);
    end
end