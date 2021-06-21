function S = cpmat(s)
    S = [0 -s(3) s(2);
         s(3) 0 -s(1);
         -s(2) s(1) 0];
end