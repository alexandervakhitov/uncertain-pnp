function xn = normalize_points(x, K)
    n = size(x, 2);
    xn = zeros(2, n);
    for i = 1:n
        xnh = K\[x(:, i);1];
        xn(:,i) = xnh(1:2)/xnh(3);
    end
end