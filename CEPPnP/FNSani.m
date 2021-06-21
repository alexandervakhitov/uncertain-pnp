%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the CEPPnP 
% algorithms described in:
%
%       Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer.
%       Leveraging Feature Uncertainty in the PnP Problem. 
%       In Proceedings of BMVC, 2014. 
%
% Copyright (C) <2014>  <Luis Ferraz, Xavier Binefa, Francesc Moreno-Noguer>
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
% Luis Ferraz, CMTech-UPF, September 2014.
% luisferrazc@gmail.com,http://cmtech.upf.edu/user/62
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ resTheta, vecTheta, valTheta] = FNS(M,dM,Cu,Vn)  
    %% Fundamental Numerical Scheme
       
    nU      = size(M,1)/2;
    nV      = size(M,2);
    W       = ones(1,nU); 
    CMu2    = zeros(nV,nV*nU);
    M22      = zeros(nV,nV*nU);
    
    tones = ones(nV,nV);

    Mt  = M';
    dMt = dM';
 
    ssCu = sqrt(squeeze([sum(sum(Cu(:,:,:)))])'); %sum(sum(Cu(:,:,j))) 
    ssCu = [ssCu; ssCu];
    ssCu = repmat(ssCu(:)',nV,1);
    dMt = dMt .* ssCu;
    
    tdMt = dMt(:,1:2:2*nU);
    idx = 1:nV*nU;
    idx = reshape(idx,nV,nU);
   
    mMt = reshape(Mt,nV,2,nU);
    
    for id = 1:nU
        dMit    = tdMt(:,id);
        CMu2(:,idx(:,id)) = dMit * dMit';
        Mti = mMt(:,:,id);
        M22(:,idx(:,id)) = Mti * Mti';
    end

    %indices for kronecker product
    ma = 1;
    na = nU;
    mb = nV;
    nb = nV;
    [ia,ib] = meshgrid(1:ma,1:mb);
    [ja,jb] = meshgrid(1:na,1:nb);

    
    for it = 1:100
        % 2.- Compute      
        Vnt = Vn';
        W  = Vnt * reshape(Vnt * CMu2,nV,nU);      
        
        W = 1./W;
        W = W./norm(W);
        
        W2    = W.^2;
        
        %hessian  

        %"for" unfolding
        tW  = W(ia,ja).*tones(ib,jb);
        tW2 = W2(ia,ja).*tones(ib,jb);
         
        tN  = M22.*tW;
        N  = sum(reshape(tN,nV,nV,nU),3);
        
        tCMuW2 = CMu2 .* tW2;
        VntM2Vn  = Vnt * reshape(Vnt * M22,nV,nU);
        tVntM2Vn =  VntM2Vn(ia,ja).*tones(ib,jb); %kron(VntM2Vn,tones);
        
        tL = tVntM2Vn .* tCMuW2;
        L  = sum(reshape(full(tL),nV,nV,nU),3);
     
        [~,s,v] = svd(N-L);
       
        rv = real(v(:,end));
       
        % 4.- If Theta = Theta0 up to sign => stop
        if sum(abs(sign(Vn) - sign(rv(:,1)))) < eps
            break;
        else 
            Vn = rv(:,1);
        end
    end
    
    resTheta = Vn;
    vecTheta = v;
    valTheta = s;
end