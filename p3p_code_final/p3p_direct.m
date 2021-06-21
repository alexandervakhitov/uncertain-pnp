function [R t] = p3p_direct( U, u0 )
%
if size(u0,2) ~= 3 
    disp('not minimal case!');
    return;
end

d12 = norm(U(:,1)-U(:,2))^2;
d13 = norm(U(:,1)-U(:,3))^2;
d23 = norm(U(:,2)-U(:,3))^2;

%call GB solver
[dep1 dep2 dep3] = GB_Solver_P3P_3Variable_2Order(u0(1,1),u0(2,1),u0(1,2),u0(2,2),u0(1,3),u0(2,3),d12, d13, d23);

%3D-3D relative pose problem
index = 1;
for i = 1:length(dep1)
    if dep1(i) > 0 && dep2(i) > 0 && dep3(i) > 0
        Ucam(:,1) = dep1(i)*[u0(:,1);1];
        Ucam(:,2) = dep2(i)*[u0(:,2);1];
        Ucam(:,3) = dep3(i)*[u0(:,3);1];
        [Rtemp ttemp] = calcampose(Ucam,U);
        R(:,:,index) = Rtemp;
        t(:,index) = ttemp;
        index = index+1;
    end
end        

%choose pose with positive depth
for i = 1:size(t,2)
    proj = R(:,:,i)*U + repmat(t(:,i),1,3);
    if min(proj(3,:)) < 0
        error(i) = inf;
        continue;
    else
        err = u0(1:2,:) - proj(1:2,:)./repmat(proj(3,:),2,1);
        error(i) = sum(sum(err.*err));
    end
end
   
%select solution
rr = find(error <= 1e-3);
if isempty(rr)
    rr = (error <= min(error));
end
R = R(:,:,rr);
t = t(:,rr);
end

function [R2,t2] = calcampose(XXc,XXw)

n= length(XXc);

X= XXw;
Y= XXc;

K= eye(n)-ones(n,n)/n;

ux= mean(X,2);
uy= mean(Y,2);
sigmx2= mean(sum((X*K).^2));

SXY= Y*K*(X')/n;
[U, D, V]= svd(SXY);
S= eye(3);
if det(SXY) < 0
    S(3,3)= -1;
end

R2= U*S*(V');
c2= trace(D*S)/sigmx2;
t2= uy-c2*R2*ux;

X= R2(:,1);
Y= R2(:,2);
Z= R2(:,3);
if norm(xcross(X,Y)-Z) > 2e-2
    R2(:,3)= -Z;
end

end

function c = xcross(a,b)

c = [a(2)*b(3)-a(3)*b(2);
     a(3)*b(1)-a(1)*b(3);
     a(1)*b(2)-a(2)*b(1)];
 
end
