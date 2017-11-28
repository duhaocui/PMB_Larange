function [r,x,P,cupd,rupd,xupd,Pupd,lupd] = dataAssocNew2(cupd,rupd,xupd,Pupd,aupd,lupd,cnew,rnew,xnew,Pnew,lnew)

Hpre = length(cupd);        % num of single target hypotheses updating pre-existing tracks
Hnew = length(cnew);        % num of single target hypotheses updating new tracks
H = Hpre+Hnew;              % total num of single target hypotheses
mcur = Hnew/2;              % num of measurements in current scan

options = optimoptions('intlinprog','Display','off');
% If there is no pre-existing track
if Hpre == 0
    A0 = zeros(mcur,H);
    % for each new track
    idx = 0;
    for i = 1:mcur
        A0(i,idx+1:idx+2) = 1;
        idx = idx+2;
    end
    b0 = ones(mcur,1);
    u = intlinprog(cnew,1:length(cnew),[],[],A0,b0,zeros(length(cnew),1),ones(length(cnew),1),[],options);
    I = u==1;
    r = rnew(I);
    x = xnew(:,I);
    P = Pnew(:,:,I);
    return;
end

% n: num of pre-existing tracks, if two single target hypotheses have the
% same length of trajectory and the same ancestor, they belong to the same
% track; num of new tracks = num of meas at current scan
npre = 1;
% number of single target hypotheses in pre-exist track i, i belongs {1,...,npre}
ns = ones(1,1);
% calculate the length of each trajectory, used for determine the num of
% scans should be used
tralen = zeros(Hpre,1);
tralen(1) = length(lupd{1});

if Hpre >= 2
    for i = 2:Hpre
        tralen(i) = length(lupd{i});
        if length(lupd{i-1}) == length(lupd{i}) && lupd{i-1}(1) == lupd{i}(1)
            ns(npre) = ns(npre)+1;
        else
            npre = npre+1;
            ns(npre) = 1;
        end
    end
end

% construct binary indicator matrix for constraint (1): each track should
% only be used once, each measurement received in current scan creates a
% new track
A0 = zeros(npre+mcur,H);
% for each pre-existing track
idx = 0;
for i = 1:npre
    A0(i,idx+1:idx+ns(i)) = 1;
    idx = idx+ns(i);
end
% for each new track
for i = 1:mcur
    A0(npre+i,idx+1:idx+2) = 1;
    idx = idx+2;
end

% construct binary indicator matrix for constraint (2): each measurement in
% each scan should only be used once
A2 = zeros(mcur,H); % current scan
tcur = zeros(H,1);  % target trajectory of the current scan
for i = 1:Hpre
    tcur(i) = lupd{i}(end);
end
for i = Hpre+1:H
    tcur(i) = lnew{i-Hpre};
end
for i = 1:mcur
    A2(i,tcur==i) = 1;
end

% constriant equality
b0 = ones(npre+mcur,1);
b2 = ones(mcur,1);

% cost of single target hypotheses
c = [cupd;cnew];

% the maximum length of the trajectory should be no less than the num of
% scans considered in the data association
if max(tralen) == 1
    % no need to do Lagrange relaxation
    u = intlinprog(c,1:length(c),[],[],[A0;A2],[b0;b2],zeros(length(c),1),ones(length(c),1),[],options);
    
    % single target hypotheses in the "best" global association hypotheses updating
    % pre-existing tracks
    I = u(1:Hpre)==1;
    r = rupd(I);
    x = xupd(:,I);
    P = Pupd(:,:,I);
    
    % single target hypotheses in the "best" global association hypotheses
    % updating new tracks
    I = u(Hpre+1:end)==1;
    r = [r;rnew(I)];
    x = [x xnew(:,I)];
    P = cat(3,P,Pnew(:,:,I));
    % no need to prune
    return;
end

% get number of measurements of last scan that are used to update single
% target hypotheses after pruning, extracted from trajectory
tpre = zeros(Hpre,1);         % target trajectory of the lastest scan
for i = 1:Hpre
    tpre(i) = lupd{i}(end-1); % assume that we start from time step 2 or after
end
% num of measurements in last scan, do not count missed detection and
% non-exist new track 
measUnique = unique(tpre(tpre~=0));
mpre = length(measUnique);

% construct binary indicator matrix for constraint (2): each measurement in
% each scan should only be used once
A1 = zeros(mpre,H); % last scan
for i = 1:mpre
    A1(i,[tpre==measUnique(i);false(Hnew,1)]) = 1;
end

% constriant equality
b1 = ones(mpre,1);

%%
% dual decomposition, solution is a binary indicator vector, decides which
% single target hypotheses are included in the "best" global hypotheses.
% options = optimoptions('intlinprog','Display','off');

% subproblem 1: min(c/2+\delta1)*u, s.t. [A0;A1]*u = [b0;b1];
% subproblem 2: min(c/2+\delta2)*u, s.t. [A0;A2]*u = [b0;b2];

% Larange multiplier \delta is initialised with 0
delta1 = zeros(H,1);
delta2 = zeros(H,1);

% initialise gap, maximum num of iteration
gap = inf;
numIteration = 0;
maxIteration = 2e2;

while (gap>0.1 && numIteration<maxIteration)
    % get suboptimal solution for each subproblem
    u1_hat = intlinprog((c/2+delta1),1:length(c),[],[],[A0;A1],[b0;b1],zeros(length(c),1),ones(length(c),1),[],options);
    u2_hat = intlinprog((c/2+delta2),1:length(c),[],[],[A0;A2],[b0;b2],zeros(length(c),1),ones(length(c),1),[],options);
    
    % calculate step size used in subgradient methods
    % first calculate subgradient
    g1 = u1_hat - (u1_hat+u2_hat)/2;
    g2 = u2_hat - (u1_hat+u2_hat)/2;
    % second calculate best primal cost
    [bestPrimalCost,idx] = min([c'*u1_hat,c'*u2_hat]);
    % third calculate dual cost
    dualCost = (c/2+delta1)'*u1_hat + (c/2+delta2)'*u2_hat;
    % fourth calculate step size used in subgradient method
    if ~all(g1)
        stepSize1 = 0;
    else
        stepSize1 = (bestPrimalCost - dualCost)/(norm(g1)^2);
    end
    if ~all(g2)
        stepSize2 = 0;
    else
        stepSize2 = (bestPrimalCost - dualCost)/(norm(g2)^2);
    end
    
    % update Lagrange multiplier
    delta1 = delta1 + stepSize1*g1;
    delta2 = delta2 + stepSize2*g2;
    % calculate duality gap
    gap = (bestPrimalCost - dualCost)/bestPrimalCost;
    
    numIteration = numIteration+1;
end

% output "best" assignment when meet termination criteria
if idx == 1
    u = u1_hat;
else
    u = u2_hat;
end

% single target hypotheses in the "best" global association hypotheses updating
% pre-existing tracks
I = u(1:Hpre)==1;
r = rupd(I);
x = xupd(:,I);
P = Pupd(:,:,I);

% get the length of single target hypotheses updating pre-existing tracks
prelen = length(r);

% Searching for ancester labels and trajectories of selected hypotheses, used for pruning
a = aupd(I);
l = lupd(I);

% single target hypotheses in the "best" global association hypotheses 
% updating new tracks
I = u(Hpre+1:end)==1;
r = [r;rnew(I)];
x = [x xnew(:,I)];
P = cat(3,P,Pnew(:,:,I));

%%
% N(one)-scan pruning, used to remove unrelated single target
% hypotheses from pre-existing tracks

% Besides, tracks are still need to be pruned. If a single target
% hypothesis selected in the global hypotheses has consecutive 0
% label, i.e., missed detection or non-exist, the whole track is pruned.
idx_delete = [];
for i = 1:prelen
    if length(l{i}) >=2 && r(i) <= 1e-4
        idx_delete = [idx_delete;i];
    end
end
a = a(idx_delete);
I = ismember(aupd,a); % get index of elements that can be retained
rupd = rupd(I);
xupd = xupd(:,I);
Pupd = Pupd(:,:,I);
lupd = lupd(I);
cupd = cupd(I);

end