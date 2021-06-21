function d2R = rodrigues2(rv,Q)
%     rv = rodrigues(R);
    a = rv(1);
    b = rv(2);
    c = rv(3);
    [R, dR, dRa, dRb, dRc, da1, da2, da3] = opnp_lines.rod2d(a,b,c);
    d1 = 2*(dR*Q*dRa(:)+da1*Q*R(:));
    d2 = 2*(dR*Q*dRb(:)+da2*Q*R(:));
    d3 = 2*(dR*Q*dRc(:)+da3*Q*R(:));
    d2R = [d1 d2 d3];
end