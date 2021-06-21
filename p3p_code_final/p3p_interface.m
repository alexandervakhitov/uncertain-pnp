function [R t] = p3p_interface( U, u0 )
%
if size(u0,2) ~= 3 
    disp('not minimal case!');
    return;
end

if size(u0,1) < 3
    u = [u0; ones(1,3)];
else
    u = u0;
end

for i = 1:3
    u(:,i) = u(:,i)/norm(u(:,i));
end

%call p3p solver
sol = p3p(U,u);

%extract pose
for i = 1:size(sol,2)/4
    R(:,:,i) = sol(:,4*(i-1)+2:4*i).';
    t(:,i) = -R(:,:,i)*sol(:,4*(i-1)+1);
end  

%choose pose with positive depth
for i = 1:size(sol,2)/4
    %removing complex solutions
    Rtemp = R(:,:,i); 
    ttemp = t(:,i);
    if sum(not(imag(Rtemp(:)))) < 9 || sum(not(imag(ttemp))) < 3
        error(i) = inf;
        continue;
    end
    
    %removing those behind the camera
    proj = Rtemp*U + repmat(ttemp,1,3);
    if min(proj(3,:)) < 0
        error(i) = inf;
        continue;
    else         
        err = u0(1:2,:) - proj(1:2,:)./repmat(proj(3,:),2,1);
        error(i) = sum(sum(err.*err));
    end
end
   
%select solution
rr = find(error <= 1e10);
R = R(:,:,rr);
t = t(:,rr);
end
