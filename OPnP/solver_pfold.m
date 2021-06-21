function [sols,stats] = solver_pfold(C0,settings)


if isfield(settings,'debug') == 0;
    settings.debug = 0;   % compute statistics conditioning etc. much slower
end
if isfield(settings,'nrun') == 0
    settings.nrun = 0; % # of run to find best sub-determinants in p-fold integer solver
end

if isfield(settings,'old_zplineq') == 0
    settings.old_zplineq = 0; % use old zp linear equation solver, ~slower
end

if isfield(settings,'cutoff_threshold') == 0;
    settings.cutoff_threshold = 1e12;
end

if isfield(settings,'full_QR1') == 0;
    settings.full_QR1 = 0;   % 0 : sparse QR on the first part of QR - faster and higher accuracy
end

if isfield(settings,'real') == 0;
    settings.real = 0;       % only output real solutions
end

if isfield(settings,'inverter_threshold') == 0
    settings.inverter_threshold = 1; % set this to be smaller seem to gain accuracy, but might give large error for cases with large scale differences
end

debug = settings.debug;

C0 = C0(:,settings.reorder);

C0 = bsxfun(@rdivide,C0,sqrt(sum(C0.^2,2)));


I = settings.I;
J = settings.J;
mon = settings.mon;
% keyboard
neq  = size(C0,1);
nmon = size(mon,2);

T = settings.T;
T0 = settings.Tid;
C1 = C0';
vv = (C1(T));
C=(sparse(settings.II,settings.JJ,vv));

% T = settings.T;
% C = zeros(max(I{end}),nmon);
% T0 = settings.Tid;
% C1 = C0';
% C(T0) = (C1(T));
% C = sparse(C);


P = settings.P;  % permissible
R = settings.R;  % reducible
E = settings.E;  % excessive

p = neq; %# of unknonws : assume # of unknown = # of eqs
n = nmon; %# of monomials in the template
r = settings.dim; %dim of p-fold solutions
ind = settings.ind; % monomial indices action_variable * permissibles

stats.neq = size(C,1);
stats.nmon = size(C,2);

%
% construct modulo matrix
%
% reorder
ne = length(E);
nr = length(R);
np = length(P);
if debug
    modstats.n_to_reduce = nr;
    modstats.n_excessive = ne;
end
V = [E R P];
C = C(:, [E, R, P]);

% eliminate excessive monomials (done by lu in the qr paper.
% this is more general)

% the row_echelon method was used before, but the qr version is almost
% always much faster.
if(~isempty(E))
    [qq rr ee] = qr(C(:, 1 : length(E)));
    %     CC = sparse(C(:, length(E) + 1 : end));
    C2 = [rr qq'*C(:, length(E) + 1 : end)];
    kk = abs(rr(1))./abs(diag(rr)) < settings.cutoff_threshold;
    k = find(diff(kk) == -1);
    if(isempty(k))
        k = length(kk);
    end
else
    C2 = C;
    k = 0;
end

% partition C into R- and P-parts
CR = C2(k + 1 : end, ne + 1 : ne + nr);
CP = C2(k + 1 : end, end - np + 1 : end);
mm = size(CR, 1);
nn = size(CR, 2) + size(CP, 2);
if(nn - mm > r)
    error('not enough equations for that solution dimension');
end

% eliminate R-monomials (this step is included in the lu factorization
% in the paper. qr is slightly slower but more stable).
[q2 UR2_0] = qr(full(CR));
CP = q2'*CP;

% [LL,UU,PP] = lu(C(:,1:(ne+nr)));


% select basis (qr + column pivoting)
CP2 = CP(1 : nr, :);
CP3 = CP(nr + 1 : end, :);
[q3 r3 e] = qr(CP3, 0);
CP4 = CP2(:, e(1 : end - r));
CB1 = CP2(:, e(end - r + 1 : end));
UP3 = r3(1 : np - r, 1 : np - r);
CB2 = r3(1 : np - r, end - r + 1 : end);
if(isempty(CP4)), CP4 = []; end;
if(isempty(UP3)), UP3 = []; end;
ee = [1 : ne + nr e + ne + nr];

if debug
    stats.basis = P(e(end-r+1:end));
end

V = V(ee);
% mon = mon(ee);

%
% elimination step
%
Celim = [UR2_0(1 : nr + np - r, :) [CP4; UP3]];


T = - Celim \ [CB1; CB2];


if debug
    modstats.rankdiff = size(Celim, 2) - rank(Celim);
    modstats.condition = cond(Celim);
end
% modulo matrix
modM = zeros(r, n);
modM(:, end - r + 1 : end) = eye(r);
modM(:, ne + 1 : end - r) = T';


e = V;

%
% save some statistics
%
if debug
    stats.n_monomials = nmon;
    stats.n_eqs = size(C, 1);
    stats.rank = rank(full(C));
    stats.rankdiff = size(C, 2) - stats.rank;
    stats.inner_rankdiff = modstats.rankdiff;
    stats.condition = modstats.condition;
    stats.n_permissible = length(P);
    %     stats.n_to_reduce = modstats.n_to_reduce;
    %     stats.n_excessive = modstats.n_excessive;
    %     stats.reducible_range = [stats.n_excessive + 1, stats.n_monomials];
    %     stats.reducibles = stats.reducible_range(1) : stats.reducible_range(2);
    stats.basis = mon(end - r + 1 : end);
end
%
% construct action matrix
%
% m = construct_actionmatrix(mon, modM, settings.dim, P, x);
ind2 = zeros(np,1);
ind2(P) = ind;

% [nouse,ind2] = ismember(ind2(e(n-r+1:n)),e);

ind2 = find_id(e,ind2(e(n-r+1:n)));

% ind2 = find_id(e,[ind2(e(n-r+1:n));settings.ids_a1;settings.ids_a3;settings.ids_abc]');
%
% ids_vv = ind2(end-12+1:end);
% ind2 = ind2(1:end-12);


M = zeros(n, r);
M(ind2, :) = eye(r);
m = modM * M;

[vv, dk] = eig(m');
d=diag(dk);


%% Extract possible solutions

%setup matrix
okmon = find(sum(modM,1)~=0);
% mon = mon(:,e);
% MM = mon(:,okmon)';

%setup vv
vv = modM(:,okmon)'*vv;

% PRid = e(end-nr-np+1:end);
% PRid = PRid(okmon-ne);


% ids_abc = find_id(e,settings.ids_abc)'- ne;

if strcmp(settings.p_inverter,'all') == 0 && strcmp(settings.p_inverter,'best') == 0
    sid  = settings.p_inverter;
    sid_r = 1:p;
    sid_r(sid)=[];
    
    %     [nouse,ids_a13] = ismember([settings.ids_a1(sid),settings.ids_a3(sid)],PRid);
    ids_a13 = find_id(e,[settings.ids_a1(sid),settings.ids_a3(sid)]) - ne;
    ids  = [ids_a13,ids_abc([sid,sid_r])'];
    
    vv1 = sqrt(vv(ids(2),:)./vv(ids(1),:));
    
    constant = vv(ids(1),:)./vv1;
    
    
    if settings.real
        realid = imag(vv1) == 0;
        vv1 = vv1(:,realid);
        constant = constant(realid);
    else
        realid = 1:settings.dim;
    end
    vv1_ = -vv1;
    mmid = settings.p_inverter;
else
    %     ids_ac2_list = find_id(e,settings.ids_ac2) - ne;
    
    ids_a13_list = find_id(e,[settings.ids_a1;settings.ids_a3]) - ne;
    
%     ids_a2c_list = find_id(e,[settings.ids_a2c]) - ne;
    
    ids_a1       = ids_a13_list(1:4);
    
    sols = [];
    for ii = 1:p
        ids_a13 = ids_a13_list([ii,p+ii]);
        
        sid  = ii;
        sid_rl{ii} = 1:p;
        sid_rl{ii}(sid)=[];
 
        % 
        ids  = [ids_a13];
        
        idsl{ii} = ids;
        
        % a3/a1
        vv1l{ii} = sqrt(vv(ids(2),:)./vv(ids(1),:));
        
        
        realid = imag(vv1l{ii}) == 0;
        vv1rl{ii} = vv1l{ii}(:,realid);
        
        
        constant{ii} = vv(ids(1),:)./vv1l{ii};
        
        
        fail_flag = isempty(vv1rl{ii});
        if ~fail_flag && min(abs(real((vv1l{ii})))) > 1e-7
            mm(ii) = min(norm(vv1l{ii}));
            mx(ii) = max(abs(vv1rl{ii}));
        else
            mm(ii) = inf;
            mx(ii) = 1;
        end
        
        
        if strcmp(settings.p_inverter,'all') == 1;
            % %             ids_a2c = ids_a2c_list(3*(ii-1)+(1:3));
            vv2l{ii}  =(vv(ids_a1(sid_rl{ii}),realid)./(ones(3,1)*(constant{ii}(:,realid))));
            solsl = [ [vv1rl{ii} -vv1rl{ii}]; [vv2l{ii} -vv2l{ii}] ] ;
            solsl([sid,sid_rl{ii}],:) = solsl;
            sols = [sols solsl];
        end
        
        
        
    end
    
    if strcmp(settings.p_inverter,'best') == 1
        [mme,mmid] = min(mm);
        [mxe,mxid] = max(mx);
        
        %     mm
        
        if mme < settings.inverter_threshold;
            [mme,mmid] = max(mm);
        else
            %         [nouse,mmid] = max(mx);
            mmid = 1;
        end
        
        %     [min(mm) mmid mm]
        
        vv1 = vv1l{mmid};
        cc  = constant{mmid};
        if settings.real
            realid   = imag(vv1) == 0;
            vv1 = vv1(realid);
            cc  = cc(realid);
            
        else
            realid = 1:length(vv1);
        end
        
        
        vv1_ = -vv1;
        
        
        sid = mmid;
        sid_r = sid_rl{mmid};
        ids = idsl{mmid};
        constant = constant{mmid}(realid);
        
        
    end
end

if ~strcmp(settings.p_inverter,'all');
    if 1
        vv2  =(vv(ids_a1(sid_rl{mmid}),realid)./(ones(3,1)*(constant)));
    end
    
    vv2_  = -vv2;
    
    sols = [ [vv1 vv1_]; [vv2 vv2_]];
    sols([sid,sid_r],:) = sols;
end

if settings.real
    realid = (sum(imag(sols)<1e-4)==4);
    sols   = sols(:,realid);
end




end % function end

function ids_reorder = find_id (list,ids)
% slightly faster than ismember
ccc = zeros(1,length(list));
ccc(ids) = 1:length(ids);
ccc      = ccc(list);
[nouse,ids_reorder,ff] = find(ccc);
ids_reorder(ff) = ids_reorder;

end
