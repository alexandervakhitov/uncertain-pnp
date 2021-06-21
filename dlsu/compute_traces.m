function [sigmas2] = compute_traces(Sigmas)
    d = size(Sigmas, 1);
    n = size(Sigmas, 3);
    sigmas2 = zeros(n, 1);
    for i = 1:n
        avg_eigenval = 1.0/d*trace(Sigmas(:,:, i));        
        sigmas2(i) = avg_eigenval;
    end
end