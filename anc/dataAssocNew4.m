function [r,x,P,cupd,rupd,xupd,Pupd,lupd,lambdau,xu,Pu,aupd] = dataAssocNew4...
    (cupd,rupd,xupd,Pupd,lupd,cnew,rnew,xnew,Pnew,lnew,lambdau,xu,Pu,aupd,it)

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
Htl = zeros(maxtralen,1);

A{1} = At;
b{1} = bt;
Htl(1) = H;

for tl = 2:maxtralen
    % get number of single target hypotheses existed at current-tl+1 scan
    Htemp = length(find(tralen>=tl));
    Htl(tl) = Htemp;
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
b = b(~cellfun('isempty',b));
maxtralen = length(A);

options = optimoptions('intlinprog','Display','off');
% u = intlinprog(c,1:length(c),[],[],[A0;At],[b0;bt],zeros(length(c),1),ones(length(c),1),[],options);
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
    
    % initialise gap, maximum num of iteration
    numIteration = 0;
    maxIteration = 2e2;
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
                    c_hat_nonegative = c_hat - min(c_hat(:));
                    % get number of measurements at scan tl
                    num_meas = size(A{tl},1);
                    % construct track to measurement assignment matrix at scan tl
                    cost = ones(npre+mcur,num_meas)*1e3;
                    % store assignment index of the single target hypothesis with the
                    % minimum cost in each track
                    idxCost = zeros(npre+mcur,num_meas);
                    for j = 1:num_meas
                        idx = 0;
                        for i = 1:npre
                            % find single target hypotheses in track i that use this
                            % measurement
                            is = find(A{tl}(j,idx+1:idx+ns(i))==1);
                            if isempty(is)
                                % Not found
                                cost(i,j) = 1e3;
                            else
                                % if found, find the single target hypothesis with the
                                % minimum cost, and record its index
                                [cost(i,j),idxmin] = min(c_hat_nonegative(idx+is));
                                idxCost(i,j) = idx+is(idxmin);
                            end
                            idx = idx+ns(i);
                        end
                        for i = npre+1:npre+mcur
                            is = find(A{tl}(j,idx+1:idx+2)==1);
                            if isempty(is)
                                cost(i,j) = 1e3;
                            else
                                [cost(i,j),idxmin] = min(c_hat_nonegative(idx+is));
                                idxCost(i,j) = idx+is(idxmin);
                            end
                            idx = idx+2;
                        end
                    end
                    % find the most likely assignment
                    costInput = -[cost 1e6*ones(npre+mcur,npre+mcur-num_meas)];
                    costInput = costInput-min(costInput(:));
                    [assignments,~] = auctionAlgorithm(costInput);
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
                                [~,idxtratl] = min(c_hat_nonegative(missornull+idx));
                                utemp(idx+missornull(idxtratl)) = 1;
                            else
                                [~,idxmin] = min(c_hat_nonegative(idx+1:idx+ns(i)));
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
%                     utest = intlinprog((c/maxtralen+delta(:,tl)),1:length(c),[],[],...
%                         [A0;A{tl}],[b0;b{tl}],zeros(length(c),1),ones(length(c),1),[],options);
%                     if ~isequal(u_hat(:,tl),utest)
%                         1;
%                     end
            end
        end
        u_hat = round(u_hat);
        
        u_hat_mean = sum(u_hat,2)/maxtralen;
        if isempty(find(u_hat_mean~=1&u_hat_mean~=0,1))
            uprimal = u_hat(:,1);
            break;
        end
        
        % second calculate dual cost
        dualCost = sum(subDualCost);
        
        % find the best feasible primal cost
        conflictVector = u_hat_mean;
        % store conflicting single target hypotheses indices
        conflictTracks = cell(npre+mcur,1);
        % total number of possible feasible solutions
        num_feasiblePrimal = 1;
        
        idx = 0;
        % first consider pre-existing tracks
        for t = 1:npre
            if ~isempty(find(conflictVector(idx+1:idx+ns(t))~=0&conflictVector(idx+1:idx+ns(t))~=1,1))
                conflictTracks{t} = (idx+1:idx+ns(t))';
            end
            num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
            idx = idx+ns(t);
        end
        
        % second consider new tracks
        for t = npre+1:npre+mcur
            if ~isempty(find(conflictVector(idx+1:idx+2)~=0&conflictVector(idx+1:idx+2)~=1,1))
                conflictTracks{t} = (idx+1:idx+2)';
            end
            num_feasiblePrimal = num_feasiblePrimal*length(conflictTracks{t});
            idx = idx+2;
        end
        
        %%%%%%%
        Amatrix = [A0;At];
        idx_selectedHypo = conflictVector==1;
        idx_usedMeas = sum(Amatrix(:,idx_selectedHypo),2)==1;
        A_uncertain = Amatrix(~idx_usedMeas,~idx_selectedHypo);
        b_uncertain = ones(size(A_uncertain,1),1);
        len_uncertain = size(A_uncertain,2);
        c_uncertain = c(~idx_selectedHypo);
        uprimal_uncertain = intlinprog(c_uncertain,1:len_uncertain,[],[],...
            sparse(A_uncertain),b_uncertain,zeros(len_uncertain,1),ones(len_uncertain,1),[],options);
        uprimal_uncertain = round(uprimal_uncertain);
        uprimal = u_hat_mean;
        uprimal(~idx_selectedHypo) = uprimal_uncertain;
        bestPrimalCosthat = c'*uprimal;
        %%%%%%%
        
        % store indices of uncertain hypotheses
%         idx_conflictHypo = false(H,1);
%         idx_conflictHypo(cell2mat(conflictTracks)) = true;
%         idx_uncertain = idx_conflictHypo;
%         idx_certain = ~idx_uncertain;
%         A_certain = [A0;At];
%         A_uncertain = A_certain(:,idx_uncertain);
%         A_certain = A_certain(:,idx_certain);
%         b_uncertain = [b0;bt] - A_certain*u_hat_mean(idx_certain);
%         len_uncertain = size(A_uncertain,2);
%         c_uncertain = c(idx_uncertain);
%         uprimal_uncertain = intlinprog(c_uncertain,1:len_uncertain,[],[],...
%             sparse(A_uncertain),b_uncertain,zeros(len_uncertain,1),ones(len_uncertain,1),[],options);
% 
%         if isempty(uprimal_uncertain)
%             % if this reconstruction is not feasible, set conflicting new
%             % tracks to null-hypothesis
%             1;
%         else
%             uprimal_uncertain = round(uprimal_uncertain);
%             uprimal = u_hat_mean;
%             uprimal(idx_uncertain) = uprimal_uncertain;
%             bestPrimalCosthat = c'*uprimal;
%         end
        
        if bestPrimalCosthat < bestPrimalCost
            bestPrimalCost = bestPrimalCosthat;
        end
        %%%%%%
        
        gap = (bestPrimalCost - dualCost)/bestPrimalCost;
        if gap < 0.05
            break;
        end
        
        % calculate step size used in subgradient methods
        % third calculate subgradient
        g = u_hat - u_hat_mean;
        
        % fourth calculate step size used in subgradient method
        stepSize = (bestPrimalCosthat - dualCost)/(norm(g)^2)/2;
        
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

% get trajecotries of ML tracks, used for pruning (do not prune new tracks)
l = lupd(I);

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

N = 3; % N-scan pruning
idx = 0;
idx_remain = [];
for i = 1:npre
    if length(l{i}) >=3 && ns(i) == 1 && isequal(l{i}(end-2:end),zeros(3,1))
        % prune null-hypothesis
    elseif length(l{i})>=N+1
        traCompared = l{i}(1:end-N+1);
        for j = idx+1:idx+ns(i)
            if isequal(lupd{j}(1:end-N+1),traCompared)
                idx_remain = cat(2,idx_remain,j);
            end
        end
    else
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