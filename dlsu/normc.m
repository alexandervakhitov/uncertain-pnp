function A_normed = normc(A)
    A2 = A.^2;
    a_norm = sqrt(sum(A2, 1));
    A_normed = A ./ repmat(a_norm, size(A, 1), 1);
    
end