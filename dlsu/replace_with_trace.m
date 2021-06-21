function [Sigmas, sigmas2] = replace_with_trace(Sigmas)
    d = size(Sigmas, 1);
    n = size(Sigmas, 3);
    sigmas2 = zeros(n, 1);
    for i = 1:n
        avg_eigenval = 1.0/d*trace(Sigmas(:,:, i));
        Sigmas(:,:,i) = eye(d)*avg_eigenval;
        sigmas2(i) = avg_eigenval;
    end
end