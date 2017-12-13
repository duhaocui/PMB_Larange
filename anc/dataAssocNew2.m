function [r,x,P,cupd,rupd,xupd,Pupd,lupd,lambdau,xu,Pu] = dataAssocNew2(cupd,rupd,xupd,Pupd,~,lupd,cnew,rnew,xnew,Pnew,lnew,lambdau,xu,Pu)

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

% n: num of pre-existing tracks, if two single target hypotheses have the
% same length of trajectory and the same ancestor, they belong to the same
% track; num of new tracks = num of meas at current scan
npre = 1;
% number of single target hypotheses in pre-exist track i, i belongs {1,...,npre}
ns = ones(1,1);
% calculate the length of each trajectory, used for determine the num of
% scans should be used
tralen = zeros(Hpre,1);
% the first track is (one of) the oldest one, thus having the longest
% trajectory
tralen(1) = length(lupd{1});

for i = 2:Hpre
    tralen(i) = length(lupd{i});
    % consider the special case of non-exist track
    if (length(lupd{i-1}) == length(lupd{i}) && (lupd{i-1}(1) == lupd{i}(1) ...
            || lupd{i-1}(1) == 0))
        ns(npre) = ns(npre)+1;
    else
        npre = npre+1;
        ns(npre) = 1;
    end
end

% cost of single target hypotheses
c = [cupd;cnew];

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
At = zeros(mcur,H); % current scan
tcur = zeros(H,1);  % target trajectory of the current scan
for i = 1:Hpre
    tcur(i) = lupd{i}(end);
end
for i = Hpre+1:H
    tcur(i) = lnew{i-Hpre};
end
for i = 1:mcur
    At(i,tcur==i) = 1;
end

% constriant equality
b0 = ones(npre+mcur,1);
bt = ones(mcur,1);

maxtralen = tralen(1); % maxtralen-scan data association

maxtralen(maxtralen>6) = 6; % sliding window

%%%%%%%%%%
maxtralen = 2;
A1 = At;
b1 = bt;

for tl = 2:maxtralen
    % get number of measurements of current-tl+1 scan that are used to update single
    % target hypotheses after pruning, extracted from trajectory
    Htemp = length(find(tralen>=tl));
    tratemp = zeros(Htemp,1);             
    for i = 1:Htemp
        tratemp(i) = lupd{i}(end-tl+1);  % target trajectory of the current-tl+1 scan
    end
    % num of measurements in current-tl+1 scan, do not count missed detection and
    % non-exist new track
    measUnique = unique(tratemp(tratemp~=0));
    mtemp = length(measUnique);
    
    Atemp = zeros(mtemp,Htemp); % current-tl+1 scan
    for i = 1:mtemp
        Atemp(i,tratemp==measUnique(i)) = 1;
    end
    Atemp = cat(2,Atemp,false(mtemp,H-Htemp));
    btemp = ones(mtemp,1);
    
    %%%%%%%%%%
    A2 = Atemp;
    b2 = btemp;
    
    At = cat(1,At,Atemp);
    bt = cat(1,bt,btemp);
end

options = optimoptions('intlinprog','Display','off');
% u = intlinprog(c,1:length(c),[],[],[A0;At],[b0;bt],zeros(length(c),1),ones(length(c),1),[],options);

%%
if 1
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
bestPrimalCost = inf;
numIteration = 0;
maxIteration = 2e2;

% keep record of unfeasible solutions
% unfeasibleComb = zeros(0,npre+mcur);

while (gap>0.05 && numIteration<maxIteration)
    % get suboptimal solution for each subproblem
    u1_hat = intlinprog((c/2+delta1),1:length(c),[],[],[A0;A1],[b0;b1],zeros(length(c),1),ones(length(c),1),[],options);
    u2_hat = intlinprog((c/2+delta2),1:length(c),[],[],[A0;A2],[b0;b2],zeros(length(c),1),ones(length(c),1),[],options);

    % first calculate best primal cost
    % find tracks with conflicting hypotheses
    conflictVector = (u1_hat+u2_hat)/maxtralen;
    % store conflicting single target hypotheses indices
    conflictTracks = cell(npre+mcur,1);
    idx = 0;
    % total number of possible feasible solutions
    num_feasiblePrimal = 1;
    
    % first consider pre-existing tracks
    for t = 1:npre
        conflictTracks{t} = find(conflictVector(idx+1:idx+ns(t))~=0&conflictVector(idx+1:idx+ns(t))~=1)+idx;
        if isempty(conflictTracks{t})
            conflictTracks{t} = find(conflictVector(idx+1:idx+ns(t))==1)+idx;
        end
        idx = idx+ns(t);
    end
    
    for t = 1:npre
%         if length(conflictTracks{t}) > 1
%             measUsedother = cellfun(@(x) x(end-1), lupd(cell2mat(conflictTracks(setdiff(1:npre,t)))));
%             measUsedself = cellfun(@(x) x(end-1), lupd(cell2mat(conflictTracks(t))));
%             measUsedself = setdiff(measUsedself,0);
%             if isempty(intersect(measUsedself,measUsedother))
%                 % do not consider non-exist single target hypothesis
%                 if ~any(lupd{conflictTracks{t}(1)})
%                     conflictTracks{t} = setdiff(conflictTracks{t},conflictTracks{t}(1));
%                 end
%             end
%         end
        num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
    end
    
    % second consider new tracks
    selectedHypothesis = cell2mat(conflictTracks(1:npre));
    measUsed = cellfun(@(x) x(end), lupd(selectedHypothesis));
    for t = npre+1:npre+mcur
        conflictTracks{t} = find(conflictVector(idx+1:idx+2)~=0&conflictVector(idx+1:idx+2)~=1)+idx;
        if isempty(conflictTracks{t})
            conflictTracks{t} = find(conflictVector(idx+1:idx+2)==1)+idx;
        end
        if isempty(find(measUsed==t-npre,1))
            conflictTracks{t} = idx+2;
        end
        idx = idx+2;
        num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
    end
    if num_feasiblePrimal==1
        % only one feasible primal solution
        uPrimal = zeros(H,1);
        uPrimal(cell2mat(conflictTracks)) = 1;
        bestPrimalCosthat = c'*uPrimal;
        if bestPrimalCosthat < bestPrimalCost
            bestPrimalCost = bestPrimalCosthat;
        end
    else
        % find the best feasible primal solution by brute force
        bestPrimalCosthat = inf;
        combs = allcombs(conflictTracks{:});
%         if ~isempty(unfeasibleComb)
%             combs = setdiff(combs,unfeasibleComb,'rows');
%         end
        [num_possiblePrimal,~] = size(combs);
        %%%%%
        allzerolen = zeros(num_possiblePrimal,1);
        for np = 1:num_possiblePrimal
            utemp = zeros(H,1);
            utemp(combs(np,:)) = 1;
            ifallzero = [A0;At]*utemp - [b0;bt];
            allzerolen(np) = length(find(ifallzero==-1));
            if ~any(ifallzero)
                tempCost = c'*utemp;
                if tempCost < bestPrimalCosthat
                    bestPrimalCosthat = tempCost;
                    uPrimal = utemp;
                end
            else
%                 unfeasibleComb = cat(1,unfeasibleComb,combs(np,:));
            end
        end
        if bestPrimalCosthat == inf
            1;
        end
        if bestPrimalCosthat < bestPrimalCost
            bestPrimalCost = bestPrimalCosthat;
        end
    end

    % second calculate dual cost
    dualCost = (c/2+delta1)'*u1_hat + (c/2+delta2)'*u2_hat;
    
    % calculate step size used in subgradient methods
    % third calculate subgradient
    g1 = u1_hat - (u1_hat+u2_hat)/2;
    g2 = u2_hat - (u1_hat+u2_hat)/2;
    
    % fourth calculate step size used in subgradient method
    if ~any(g1)
        stepSize1 = 0;
    else
        stepSize1 = (bestPrimalCost - dualCost)/(norm(g1)^2);
    end
    if ~any(g2)
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
u = uPrimal;
end
%%

% single target hypotheses in the ML global association hypotheses updating
% pre-existing tracks
I = u(1:Hpre)==1;
r = rupd(I);
x = xupd(:,I);
P = Pupd(:,:,I);

% get trajecotries of ML tracks, used for pruning
l = lupd(I);

% single target hypotheses in the ML global association hypotheses
% updating new tracks
I = u(Hpre+1:end)==1;
r = [r;rnew(I)];
x = [x xnew(:,I)];
P = cat(3,P,Pnew(:,:,I));

% N-scan pruning, used to remove unrelated single target
% hypotheses from pre-existing tracks

% Besides, tracks are still need to be pruned. If a single target
% hypothesis selected in the ML global hypotheses has consecutive 0
% label, i.e., missed detection or non-exist, the whole track is pruned.

N = 1; % N-scan pruning (currently only works for N==1, debug???)
if min(tralen)>N
    % select out the last but N elements of each trajectory
    lupdselect = cellfun(@(x) x(end-N), lupd);
    lselect = cellfun(@(x) x(end-N), l);
else
    lupdselect = ones(Hpre,1);
    lselect = ones(npre,1);
end

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
    else
        % single target hypotheses left after N-scan pruning
        if lselect(i) > 0
            idx_remain = cat(1,idx_remain,find(lupdselect(idx+1:idx+ns(i))==lselect(i))+idx);
        else
            if l{i}(1) > 0  % missed detection
                if ~any(lupd{idx+1})
                    idx_remain = cat(1,idx_remain,setdiff(find(lupdselect(idx+1:idx+ns(i))==lselect(i)),1)+idx);
                else
                    idx_remain = cat(1,idx_remain,find(lupdselect(idx+1:idx+ns(i))==lselect(i))+idx);
                end
            else            % non-exist single target hypothesis
                idx_remain = cat(1,idx_remain,1+idx);
            end
        end
    end
    idx = idx+ns(i);
end
rupd = rupd(idx_remain);
xupd = xupd(:,idx_remain);
Pupd = Pupd(:,:,idx_remain);
lupd = lupd(idx_remain);
cupd = cupd(idx_remain);

end