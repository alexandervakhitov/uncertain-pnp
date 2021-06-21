function [xopt, fmin] = refine_pnpl_levmar(Rr0, M0)
    imvect = 0;
    x0 = rodrigues(Rr0);
    
    options=[1E-01, 1E-10, 1E-10, 1E-10, 1e-6];
    itmax = 200;
    
    [ret, popt, info, covar]=levmar('opnp_lines.pnpl_res', x0, imvect, itmax, options, M0);        
    xopt = popt;
    fmin = info(2);
end