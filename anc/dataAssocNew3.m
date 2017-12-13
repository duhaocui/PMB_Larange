function [r,x,P,cupd,rupd,xupd,Pupd,lupd,lambdau,xu,Pu,aupd] = dataAssocNew3...
    (cupd,rupd,xupd,Pupd,lupd,cnew,rnew,xnew,Pnew,lnew,lambdau,xu,Pu,aupd)

Hpre = length(cupd);        % num of single target hypotheses updating pre-existing tracks
Hnew = length(cnew);        % num of single target hypotheses updating new tracks
H = Hpre+Hnew;              % total num of single target hypotheses
mcur = Hnew/2;              % num of measurements in current scan

% If there is no pre-existing track
if Hpre == 0
    r = rnew(2:2:end);
    x = xnew(:,2:2:end);
    P = Pnew(:,:,2:2:end);
    return;
end

% n: num of pre-existing tracks; num of new tracks = num of meas at current scan
unique_a = unique(aupd,'stable');
npre = length(unique_a);
% number of single target hypotheses in pre-exist track i, i belongs {1,...,npre}
ns = zeros(npre,1);
for i = 1:npre
    ns(i) = length(find(aupd==unique_a(i)));
end

% calculate the length of trajectory of each single target hypothesis in 
% pre-existing tracks, used for determine the num of scans should be used
tralen = cellfun(@(x) length(x),lupd);

% cost of single target hypotheses
c = [cupd;cnew];
% single target hypotheses are selected, one from each track
idx = 0;
for i = 1:npre
    c(idx+1:idx+ns(i)) = c(idx+1:idx+ns(i)) - min(c(idx+1:idx+ns(i)));
    idx = idx+ns(i);
end

% construct binary indicator matrix for constraint (1): each track should
% only be used once, each measurement received in current scan creates a
% new track; this constraint is neccessary for the implementation of dual
% decomposition, since sliding window is used.
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
% target trajectory of the current scan
tcur = cat(1,cellfun(@(x) x(end),lupd),cellfun(@(x) x(1),lnew));
At = zeros(mcur,H);
for i = 1:mcur
    At(i,tcur==i) = 1;
end

% constriant equality
b0 = ones(npre+mcur,1);
bt = ones(mcur,1);

maxtralen = max(tralen); % maxtralen-scan data association

slideWindow = 4; % apply sliding window
maxtralen(maxtralen>slideWindow) = slideWindow; 

A = cell(maxtralen,1);
b = cell(maxtralen,1);

A{1} = At;
b{1} = bt;

for tl = 2:maxtralen
    % get number of single target hypotheses existed at current-tl+1 scan
    Htemp = length(find(tralen>=tl));
    % target trajectory of the current-tl+1 scan, tracks from left to
    % right, old to new
    tratemp = cellfun(@(x) x(end-tl+1),lupd(1:Htemp));
    % num of measurements in current-tl+1 scan, do not count missed detection and
    % non-exist new track
    measUnique = unique(tratemp(tratemp~=0),'stable');
    mtemp = length(measUnique);
    
    Atemp = zeros(mtemp,Htemp); % current-tl+1 scan
    for i = 1:mtemp
        Atemp(i,tratemp==measUnique(i)) = 1;
    end
    Atemp = cat(2,Atemp,zeros(mtemp,H-Htemp));
    btemp = ones(mtemp,1);
    
    A{tl} = Atemp;
    b{tl} = btemp;
    
    At = cat(1,At,Atemp);
    bt = cat(1,bt,btemp);
end

A = A(~cellfun('isempty',A));
maxtralen = length(A);

% options = optimoptions('intlinprog','Display','off');
% u = intlinprog(c,1:length(c),[],[],[A0;At],[b0;bt],zeros(length(c),1),ones(length(c),1),[],options);

%%

% dual decomposition, solution is a binary indicator vector, decides which
% single target hypotheses are included in the "best" global hypotheses.

% subproblem t: min(c/t+\deltat)*u, s.t. [A0;At]*u = [b0;bt];

% Larange multiplier \delta is initialised with 0
delta = zeros(H,maxtralen);
% subproblem solutions
u_hat = zeros(H,maxtralen);

% initialise gap, maximum num of iteration
dualCost = 0;
numIteration = 0;
maxIteration = 1e3;

while (numIteration<maxIteration)
    % get suboptimal solution for each subproblem
    for tl = 1:maxtralen
%         u_hat(:,tl) = intlinprog((c/maxtralen+delta(:,tl)),1:length(c),[],[],...
%             [A0;A{tl}],[b0;b{tl}],zeros(length(c),1),ones(length(c),1),[],options);
        c_hat = c/maxtralen+delta(:,tl);
        % get number of measurements at scan tl
        num_meas = size(A{tl},1);
        % construct track to measurement assignment matrix at scan tl 
        cost = zeros(npre+mcur,num_meas);
        % store assignment index of the single target hypothesis with the
        % minimum cost in each track
        idxCost = zeros(npre+mcur,num_meas);
        for j = 1:num_meas
            idx = 0;
            for i = 1:npre
                is = find(A{tl}(j,idx+1:idx+ns(i))==1);
                if isempty(is)
                    % Not found
                    cost(i,j) = inf;
                else
                    % if found, find the single target hypothesis with the
                    % minimum cost, and record its index
                    [cost(i,j),idxmin] = min(c_hat(idx+is));
                    idxCost(i,j) = idx+is(idxmin);
                end
                idx = idx+ns(i);
            end
            for i = npre+1:npre+mcur
                is = find(A{tl}(j,idx+1:idx+2)==1);
                if isempty(is)
                    cost(i,j) = inf;
                else
                    [cost(i,j),idxmin] = min(c_hat(idx+is));
                    idxCost(i,j) = idx+is(idxmin);
                end
                idx = idx+2;
            end
        end
        % find the most likely assignment
        [assignments,~] = murty_custom(cost-min(cost(:)),1);
        indicator = [];
        for i = 1:npre+mcur
            if assignments(i) > 0
                indicator = cat(1,indicator,idxCost(i,assignments(i)));
            end
        end
        % if a track has no measurement assigned to it, choose the single
        % target hypotheses with the minimum cost (either, non-exist or miss)
        utemp = zeros(H,1);
        utemp(indicator) = 1;
        idx = 0;
        for i = 1:npre
            if ~any(utemp(idx+1:idx+ns(i)))
                [~,idxmin] = min(c(idx+1:idx+ns(i)));
                utemp(idx+idxmin) = 1;
            end
            idx = idx+ns(i);
        end
        for i = 1:mcur
            if ~any(utemp(idx+1:idx+2))
                utemp(idx+1) = 1;
            end
            idx = idx+2;
        end
        u_hat(:,tl) = utemp;
    end
    
    u_hat_mean = sum(u_hat,2)/maxtralen;
    if isempty(find(u_hat_mean~=1&u_hat_mean~=0,1))
        break;
    end

    % second calculate dual cost
    dualCosttemp = 0;
    for tl = 1:maxtralen
        dualCosttemp = dualCosttemp + (c/maxtralen+delta(:,tl))'*u_hat(:,tl);
    end

    if dualCosttemp > dualCost
        dualCost = dualCosttemp;
    else
        break;
    end
    
    % calculate step size used in subgradient methods
    % third calculate subgradient
    g = u_hat - u_hat_mean;
    
    % fourth calculate step size used in subgradient method
    stepSize = zeros(1,maxtralen);
    for tl = 1:maxtralen
        if ~any(g(:,tl))
            stepSize(tl) = 0;
        else
            stepSize(tl) = norm(g(:,tl))/sqrt(numIteration+1);
        end
    end
    
    % update Lagrange multiplier
    delta = delta + stepSize.*g;

    numIteration = numIteration+1;
end

%%
% output most likely assignment when meet termination criteria
if isempty(find(u_hat_mean~=1&u_hat_mean~=0,1))
    u = u_hat(:,1);
else
    conflictVector = u_hat_mean;
    % store conflicting single target hypotheses indices
    conflictTracks = cell(npre+mcur,1);
    % total number of possible feasible solutions
    num_feasiblePrimal = 1;
    
    idx = 0;
    % first consider pre-existing tracks
    for t = 1:npre
        conflictTracks{t} = find(conflictVector(idx+1:idx+ns(t))~=0&conflictVector(idx+1:idx+ns(t))~=1)+idx;
        if isempty(conflictTracks{t})
            conflictTracks{t} = find(conflictVector(idx+1:idx+ns(t))==1,1)+idx;
        end
        num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
        idx = idx+ns(t);
    end
    
    selectedHypothesis = cell2mat(conflictTracks(1:npre));
    measUsed = cellfun(@(x) x(end), lupd(selectedHypothesis));
    % second consider new tracks
    for t = npre+1:npre+mcur
        conflictTracks{t} = find(conflictVector(idx+1:idx+2)~=0&conflictVector(idx+1:idx+2)~=1)+idx;
        if isempty(conflictTracks{t})
            conflictTracks{t} = find(conflictVector(idx+1:idx+2)==1,1)+idx;
        end
        % use the constraint that measurements at current scan must be used
        % to remove unfeasible primal solution
        if isempty(find(measUsed==t-npre,1))
            conflictTracks{t} = idx+2;
        end
        num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
        idx = idx+2;
    end
    
    % If there are too many conflicts, use branch & bound instead
%     if num_feasiblePrimal < 1e4
        bestPrimalCosthat = inf;
        u = u_hat(:,1);
        combs = allcombs(conflictTracks{:});
        for np = 1:num_feasiblePrimal
            utemp = zeros(H,1);
            utemp(combs(np,:)) = 1;
            ifallzero = [A0;At]*utemp - [b0;bt];
            if ~any(ifallzero)
                tempCost = c'*utemp;
                if tempCost < bestPrimalCosthat
                    bestPrimalCosthat = tempCost;
                    u = utemp;
                end
            end
        end
%     else
%         options = optimoptions('intlinprog','Display','off');
%         u = intlinprog(c,1:length(c),[],[],[A0;At],[b0;bt],zeros(length(c),1),ones(length(c),1),[],options);
%         if isempty(u)
%             u = u_hat(:,1);
%         end
%     end
    
end

%%

% single target hypotheses in the ML global association hypotheses updating
% pre-existing tracks
I = u(1:Hpre)==1;
r = rupd(I);
x = xupd(:,I);
P = Pupd(:,:,I);

% get trajecotries of ML tracks, used for pruning (do not prune new tracks)
l = lupd(I);

% single target hypotheses in the ML global association hypotheses
% updating new tracks
I = u(Hpre+1:end)==1;
r = [r;rnew(I)];
x = [x xnew(:,I)];
P = cat(3,P,Pnew(:,:,I));

% N-scan pruning, used to remove unrelated single target
% hypotheses from pre-existing tracks. In other words, for each track, only
% single target hypotheses that have the same trajectory from time step 1
% to current - N + 1 remain

% Besides, tracks are still need to be pruned. If a single target
% hypothesis selected in the ML global hypotheses has consecutive 0
% label, i.e., missed detection or non-exist, the whole track is pruned.

N = 3; % N-scan pruning
idx = 0;
idx_remain = [];
for i = 1:npre
    if length(l{i})>=2 && r(i) == 0
        % tracks to be pruned
    elseif length(l{i})>=2 && r(i)<1e-1
        % tracks to be recycled (pruned from MBM)
        lambdau = cat(2,lambdau,r(i)');
        xu = cat(2,xu,x(:,i));
        Pu = cat(3,Pu,P(:,:,i));
    elseif length(l{i})>=N+1
        traCompared = l{i}(1:end-N);
        for j = idx+1:idx+ns(i)
            if lupd{j}(1:end-N) == traCompared
                idx_remain = cat(2,idx_remain,j);
            end
        end
    else
        % all kept
        idx_remain = cat(2,idx_remain,idx+1:idx+ns(i));
    end
    idx = idx+ns(i);
end

rupd = rupd(idx_remain);
xupd = xupd(:,idx_remain);
Pupd = Pupd(:,:,idx_remain);
lupd = lupd(idx_remain);
cupd = cupd(idx_remain);
aupd = aupd(idx_remain);

end