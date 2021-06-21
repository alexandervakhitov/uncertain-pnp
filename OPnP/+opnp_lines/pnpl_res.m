function res = pnpl_res(x, M)
    R = rodrigues(x);
    r = R(:);
    res = r'*M*r;
end