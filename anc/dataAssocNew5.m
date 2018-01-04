function [r,x,P,cupd,rupd,xupd,Pupd,lupd,lambdau,xu,Pu,aupd,cnew,rnew,xnew,Pnew,lnew,anew] = dataAssocNew5...
    (cupd,rupd,xupd,Pupd,lupd,cnew,rnew,xnew,Pnew,lnew,anew,lambdau,xu,Pu,aupd)

Hpre = length(cupd);        % num of single target hypotheses updating pre-existing tracks
Hnew = length(cnew);        % num of single target hypotheses updating new tracks
H = Hpre+Hnew;              % total num of single target hypotheses
mcur = Hnew/2;              % num of measurements in current scan

% If there is no pre-existing track
if Hpre == 0
    r = rnew;
    x = xnew;
    P = Pnew;
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

slideWindow = 6; % apply sliding window
maxtralen(maxtralen>slideWindow) = slideWindow;

A = cell(maxtralen,1);
b = cell(maxtralen,1);
% Htl = zeros(maxtralen,1);

A{1} = At;
b{1} = bt;
% Htl(1) = H;

for tl = 2:maxtralen
    % get number of single target hypotheses existed at current-tl+1 scan
    idx = find(tralen>=tl);
    Htemp = length(idx);
%     Htl(tl) = Htemp;
    % target trajectory of the current-tl+1 scan, tracks from left to
    % right, old to new
    tratemp = cellfun(@(x) x(end-tl+1),lupd(idx));
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
b = b(~cellfun('isempty',b));
maxtralen = length(A);

Amatrix = [A0;At];

options = optimoptions('intlinprog','Display','off');
% u = intlinprog(c,1:length(c),[],[],Amatrix,[b0;bt],zeros(length(c),1),ones(length(c),1),[],options);
% u = round(u);

%%
if 1
    % dual decomposition, solution is a binary indicator vector, decides which
    % single target hypotheses are included in the "best" global hypotheses.
    
    % subproblem t: min(c/t+\deltat)*u, s.t. [A0;At]*u = [b0;bt];
    
    % Larange multiplier \delta is initialised with 0
    delta = zeros(H,maxtralen);
    % subproblem solutions
    u_hat = zeros(H,maxtralen);
    
    % initialise maximum num of iteration
    numIteration = 0;
    maxIteration = 1e2;
    % store the best feasible primal cost obtained so far (upper bound)
    bestPrimalCost = inf;
    
    while (numIteration<maxIteration)
        % get suboptimal solution for each subproblem
        subDualCost = zeros(maxtralen,1);
        for tl = 1:maxtralen
            method = 2;
            switch method
                case 1
                    % implementation using branch and bound
                    [u_hat(:,tl),subDualCost(tl)] = intlinprog((c/maxtralen+delta(:,tl)),1:length(c),[],[],...
                        [A0;A{tl}],[b0;b{tl}],zeros(length(c),1),ones(length(c),1),[],options);
                case 2
                    % implementation using auction
                    c_hat = c/maxtralen+delta(:,tl);
                    % get number of measurements at scan tl
                    num_meas = size(A{tl},1);
                    % construct track to measurement assignment matrix at scan tl
                    cost = ones(npre+mcur,num_meas)*1e4;
                    % store assignment index of the single target hypothesis with the
                    % minimum cost in each track
                    idxCost = zeros(npre+mcur,num_meas);
                    for j = 1:num_meas
                        idx = 0;
                        for i = 1:npre
                            % find single target hypotheses in track i that use this
                            % measurement
                            is = find(A{tl}(j,idx+1:idx+ns(i))==1);
                            if ~isempty(is)
                                % if found, find the single target hypothesis with the
                                % minimum cost, and record its index
                                [cost(i,j),idxmin] = min(c_hat(idx+is));
                                idxCost(i,j) = idx+is(idxmin);
                            end
                            idx = idx+ns(i);
                        end
                        for i = npre+1:npre+mcur
                            is = find(A{tl}(j,idx+1:idx+2)==1);
                            if ~isempty(is)
                                [cost(i,j),idxmin] = min(c_hat(idx+is));
                                idxCost(i,j) = idx+is(idxmin);
                            end
                            idx = idx+2;
                        end
                    end
                    % find the most likely assignment
                    costInput = -[cost 1e4*ones(npre+mcur,npre+mcur-num_meas)];
                    costInput = costInput-min(costInput(:));
                    [assignments,~] = auctionAlgorithm_mex(costInput);
                    assignments(assignments>num_meas) = 0;
                    indicator = [];
                    for i = 1:npre+mcur
                        if assignments(i) > 0
                            indicator = cat(1,indicator,idxCost(i,assignments(i)));
                        end
                    end
                    % if a track has no measurement assigned to it, choose the single
                    % target hypotheses to be non-exist or miss if the track exists
                    % before scan N-tl, (if exists after scan N-tl, ofcourse, no measurement assigned)
                    utemp = zeros(H,1);
                    utemp(indicator) = 1;
                    idx = 0;
                    for i = 1:npre
                        if ~any(utemp(idx+1:idx+ns(i)))
                            if length(lupd{idx+1}) >= tl
                                tratl = cellfun(@(x) x(end-tl+1),lupd(idx+1:idx+ns(i)));
                                missornull = find(tratl==0);
                                [~,idxtratl] = min(c_hat(missornull+idx));
                                utemp(idx+missornull(idxtratl)) = 1;
                            else
                                [~,idxmin] = min(c_hat(idx+1:idx+ns(i)));
                                utemp(idx+idxmin) = 1;
                            end
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
                    subDualCost(tl) = c_hat'*u_hat(:,tl);
            end
        end
        u_hat = round(u_hat);
        
        u_hat_mean = sum(u_hat,2)/maxtralen;
        % All the subproblem solutions are equal means we have found the
        % optimal solution
        if isempty(find(u_hat_mean~=1&u_hat_mean~=0,1))
            uprimal = u_hat(:,1);
            break;
        end
        
        % second calculate dual cost
        dualCost = sum(subDualCost);
        
        % find a feasible primal solution using branch&bound
        idx_selectedHypo = u_hat_mean==1;
        idx_unselectedHypo = ~idx_selectedHypo;
        idx_unusedMeas = sum(Amatrix(:,idx_selectedHypo),2)==0;
        % if we are certain about the track, remove the single target
        % hypotheses of it
        idx = 0;
        for i = 1:npre
            if idx_unusedMeas(i) == false
                idx_unselectedHypo(idx+1:idx+ns(i)) = false;
            end
            idx = idx + ns(i);
        end
        A_uncertain = Amatrix(idx_unusedMeas,idx_unselectedHypo);
        b_uncertain = ones(size(A_uncertain,1),1);
        len_uncertain = size(A_uncertain,2);
        c_uncertain = c(idx_unselectedHypo);
        uprimal_uncertain = intlinprog(c_uncertain,1:len_uncertain,[],[],...
            sparse(A_uncertain),b_uncertain,zeros(len_uncertain,1),ones(len_uncertain,1),[],options);
        uprimal = u_hat_mean;
        uprimal(idx_unselectedHypo) = round(uprimal_uncertain);
        bestPrimalCosthat = c'*uprimal;
        
        if bestPrimalCosthat < bestPrimalCost
            bestPrimalCost = bestPrimalCosthat;
        end
        
        gap = (bestPrimalCost - dualCost)/bestPrimalCost;
        if gap < 0.05
            break;
        end
        
        % calculate step size used in subgradient methods
        % third calculate subgradient
        g = u_hat - u_hat_mean;
        
        % fourth calculate step size used in subgradient method
        stepSize = (bestPrimalCosthat - dualCost)/(norm(g)^2);
        
        % update Lagrange multiplier
        delta = delta + stepSize*g;
        
        numIteration = numIteration+1;
    end
    
    numIteration
    u = uprimal;
    
end
%%

% single target hypotheses in the ML global association hypotheses updating
% pre-existing tracks
I = u(1:Hpre)==1;
r = rupd(I);
x = xupd(:,I);
P = Pupd(:,:,I);
a = aupd(I);

% get trajecotries of ML tracks, used for pruning (do not prune new tracks)
l = lupd(I);

% recycling
% idx = r~=0&r<0.05;
% lambdau = [lambdau r(idx)'];
% xu = [xu x(:,idx)];
% Pu = cat(3,Pu,P(:,:,idx));
% r = r(~idx);
% x = x(:,~idx);
% P = P(:,:,~idx);
% a = a(idx);

% single target hypotheses in the ML global association hypotheses
% updating new tracks
I = u(Hpre+1:H)==1;
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

N = 4; % N-scan pruning
idx = 0;
idx_remain = [];
nc = 2; % prune null-hypothesis with length no less than nc
for i = 1:npre
    if (length(l{i})>=nc && ns(i)==1 && isequal(l{i}(end-nc+1:end),zeros(nc,1)))
        % prune null-hypothesis
    else
        if length(l{i})>=N+1
            traCompared = l{i}(1:end-N+1);
            for j = idx+1:idx+ns(i)
                if isequal(lupd{j}(1:end-N+1),traCompared)
                    idx_remain = cat(2,idx_remain,j);
                end
            end
        else
            idx_remain = cat(2,idx_remain,idx+1:idx+ns(i));
        end
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