function Km = compute_M_ver2_lines(U, Alph,A, Alph_S, Alph_E, l2d, wl)

% Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
% Francesc Moreno-Noguer, CVLab-EPFL, September 2007.
% fmorenoguer@gmail.com, http://cvlab.epfl.ch/~fmoreno/ 

n=size(Alph, 2); %number of 3d points
nl = size(Alph_S, 2);

fu = A(1,1);
fv = A(2,2);
u0 = A(1,3);
v0 = A(2,3);

nrows_M = 2 * (n+nl);
nrows_M = 2 * (n);
ncols_M = 12;
M=zeros(nrows_M,ncols_M);

for i=1:n
    a1=Alph(1, i);
    a2=Alph(2, i);
    a3=Alph(3, i);
    a4=Alph(4, i);
    
    ui=U(1,i);
    vi=U(2,i);
    
    %generate submatrix M
    M_=[a1*fu, 0, a1*(u0-ui), a2*fu, 0, a2*(u0-ui), a3*fu, 0, a3*(u0-ui), a4*fu, 0, a4*(u0-ui);
        0, a1*fv, a1*(v0-vi), 0, a2*fv, a2*(v0-vi), 0, a3*fv, a3*(v0-vi), 0, a4*fv, a4*(v0-vi)];
    
    
    %put M_ in the whole matrix
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=M_;
     
end

%from here works only for u0 = v0 = 0, fu = fv
% for i = 1 : nl
% 
%     as1 = Alph_S(1, i);
%     as2 = Alph_S(2, i);
%     as3 = Alph_S(3, i);
%     as4 = Alph_S(4, i);
%     lineCoefs = l2d(:, i)';
% 
%     ae1 = Alph_E(1, i);
%     ae2 = Alph_E(2, i);
%     ae3 = Alph_E(3, i);
%     ae4 = Alph_E(4, i);
% %     if (length(wl) > 1)
% %         if (length(wl) == size(l2d, 2))
% %             w = wl(i);
% %         else
% %             w = diag([wl(2*i-1), wl(2*i)]);
% %         end
% %     else
% %         w = wl;
% %     end
% %     
% %     
%     M_ = [lineCoefs*as1, lineCoefs*as2, lineCoefs*as3, lineCoefs*as4;
%           lineCoefs*ae1, lineCoefs*ae2, lineCoefs*ae3, lineCoefs*ae4];
% %       
% %     row_ini = 2*n + i*2-1;
% %     row_end = 2*n + i*2;
%     row_ini = i*2-1;
%     row_end = i*2;
%         
%     M(row_ini:row_end,:)=M_;
%     
% %     M(row_ini, 1:3) = lineCoefs*as1;
% %     M(row_ini, 4:6) = lineCoefs*as2;
% %     M(row_ini, 7:9) = lineCoefs*as3;
% %     M(row_ini, 10:12) = lineCoefs*as4;
% %     
% %     M(row_end, 1:3) = lineCoefs*ae1;
% %     M(row_end, 4:6) = lineCoefs*ae2;
% %     M(row_end, 7:9) = lineCoefs*ae3;
% %     M(row_end, 10:12) = lineCoefs*ae4;
% end


MtMe = zeros(12, 12, nl);
for i = 1:nl
    l = l2d(:, i);
    as = Alph_S(:, i);
    As = as*as';
    ae = Alph_E(:, i);
    Ae = ae*ae';
%     Mtot = zeros(12, 12);
%     for j = 1:3
%         for k = j:3
%             AddMat = l(j)*l(k)*(As + Ae);
% %             MtMe(4*(j-1)+1:4*j, 4*(k-1)+1:4*k) = MtMe(4*(j-1)+1:4*j, 4*(k-1)+1:4*k) + AddMat;            
%             Mtot(4*(j-1)+1:4*j, 4*(k-1)+1:4*k) = AddMat;
%         end
%     end
    Lt = l*l';
    MtMe(:, :, i) = kron(Lt, As + Ae);
%     Me(2*i-1, 1:4) = l(1)*as;
%     Me(2*i-1, 5:8) = l(2)*as;
%     Me(2*i-1, 9:12) = l(3)*as;
%     Me(2*i, 1:4) = l(1)*ae;
%     Me(2*i, 5:8) = l(2)*ae;
%     Me(2*i, 9:12) = l(3)*ae;
end

% for j = 1:3
%     for k = j+1:3
%         MtMe(4*(k-1)+1:4*k, 4*(j-1)+1:4*j) = MtMe(4*(j-1)+1:4*j, 4*(k-1)+1:4*k);
%     end
% end
MtMe = sum(MtMe, 3);
Mpts = M'*M;

inds = [1 5 9 2 6 10 3 7 11 4 8 12];
indsRev = [1 4 7 10 2 5 8 11 3 6 9 12];
Mpts = Mpts(indsRev, indsRev);
MtMe = MtMe + Mpts;
% M = M'*M;
%     
[V1,S1]=eig(MtMe);

Km = V1(inds,4:-1:1);

% [V,S]=eig(M);
% 
% Km = V(:,4:-1:1);
end

