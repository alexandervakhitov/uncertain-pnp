function [th, fmin] = refine_pnpl(R0, Q0, itNum)    
    [th] = opnp_lines.rodrigues(R0);    
%     itNum = 5;
        
%     errs = [];
    for i = 1:itNum
        [Ri, J] = opnp_lines.rodrigues(th);
%         H = opnp_lines.rodrigues2(th, Q0);
        r = Ri(:);
%         th = th - (J'*J + 0.00001*eye(3))\(J'*Q0*r);

%         th = th - ((J'*Q0*r)'*J'*Q0*r)\(J'*Q0*r);
%         J'
%         Q0
%         r
%         J'*Q0
%         J'*Q0*r
%         J'*Q0*r*(J'*Q0*r)'
%         J'*Q0*r
        th = th - (J'*Q0*r*(J'*Q0*r)'+0.00001*eye(3))\(J'*Q0*r);
%         th = th - 0.01*(J'*Q0*r);
%         errs = [errs; minFun(Q0, th)];
    end
    Rf = opnp_lines.rodrigues(th);
    r = Rf(:);
    fmin = r'*Q0*r;
end
function fval = minFun(Q, rv)
    R = rodrigues(rv);
    r = R(:);
    fval = r'*Q*r;
end