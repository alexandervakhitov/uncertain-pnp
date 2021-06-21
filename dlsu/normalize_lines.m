function ln = normalize_lines(l, K)
    %l'Kx_n = 0 -> ln = K'l;
    ln = K'*l;
    line_norms = sqrt(sum(ln(1:2,:).^2, 1));
    ln = ln ./ repmat(line_norms, 3, 1);
end