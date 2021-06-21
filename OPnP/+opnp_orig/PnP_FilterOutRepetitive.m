function [xx yy zz tt] = PnP_FilterOutRepetitive(xx,yy,zz,tt)
for i = 1:length(xx)
    index = [];
    for j = i+1:length(xx)
        tempi = [xx(i) yy(i) zz(i) tt(i)];
        tempj = [xx(j) yy(j) zz(j) tt(j)];
        if norm(tempi-tempj)<1e-4 || norm(tempi+tempj)<1e-4
            index = [index,j];
        end
    end
    xx(index) = [];
     yy(index) = [];
      zz(index) = [];
       tt(index) = [];
end
end