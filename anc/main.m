clc;clear
dbstop if error

% Set simulation parameters
Pd = 0.75; % probability of detection
lfai = 10; % expected number of false alarms per scan
numtruth = 6; % number of targets
simcasenum = 1; % simulation case 1 or 2 (see paper)
%simcasenum = 2;
if (simcasenum == 1) % covariance used for mid-point initialisation
    Pmid = 1e-6*eye(4);
else
    Pmid = 0.25*eye(4);
end

% Generate truth
[model,measlog,xlog] = gentruth(Pd,lfai,numtruth,Pmid,simcasenum);
load('testData7510.mat');

% Initialise filter parameters
stateDimensions = size(model.xb,1);
% Multi-Bernoulli representation
n = 0;
r = zeros(n,1);
x = zeros(stateDimensions,n);
P = zeros(stateDimensions,stateDimensions,n);
l = cell(n,1);  % store trajectory
c = zeros(n,1); % store the cost of each single target hypothesis
a = zeros(n,1);
% Unknown target PPP parameters
lambdau = model.lambdau;
xu = model.xb;
Pu = model.Pb;

% Loop through time
numTime = length(measlog);
xest = cell(numTime,1);
% GOSPA metric
gospa_vals = zeros(numTime,4);
gospa_c = 20;
gospa_p = 1;
gospa_alpha = 2;


% Two scan
for t = 1:numTime
    %%
    % Predict
%     [r,x,P,lambdau,xu,Pu] = predictStep(r,x,P,lambdau,xu,Pu,model);
    
    % Predict all single target hypotheses of previous scan
    % No need to update trajectory label, cost
    [r,x,P,lambdau,xu,Pu] = predictStep(r,x,P,lambdau,xu,Pu,model);
    
    %%
    % Update
%     [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,lupd,aupd,cupd,wnew,rnew,xnew,Pnew,lnew,cnew] = ...
%         updateStepNew(lambdau,xu,Pu,r,x,P,lpre,cpre,measlog{t},model);
    
    % Update all predicted single target hypotheses of previous scan
    [lambdau,xu,Pu,rupd,xupd,Pupd,lupd,cupd,rnew,xnew,Pnew,lnew,cnew,aupd,anew] = ...
        updateStepNew2(lambdau,xu,Pu,r,x,P,l,c,measlog{t},a,model);
    
    %%
    % Data assciation
%     [r,x,P,aselect] = dataAssocNew(wupd,rupd,xupd,Pupd,aupd,wnew,rnew,xnew,Pnew);
    
    % multi(two)-scan data association
    [r_hat,x_hat,P_hat,cupd,rupd,xupd,Pupd,lupd,lambdau,xu,Pu,aupd,cnew,rnew,xnew,Pnew,lnew,anew] = ...
        dataAssocNew5(cupd,rupd,xupd,Pupd,lupd,cnew,rnew,xnew,Pnew,lnew,anew,lambdau,xu,Pu,aupd);
    
    % Store single target hypotheses of current scan, send to next scan
    % Order: missed detection, single target hypotheses updated by the
    % first, second, third,..., measurement, followed by single target
    % hypothese of newly detected targets
    r = [rupd;rnew];
    x = [xupd xnew];
    P = cat(3,Pupd,Pnew);
    l = cat(1,lupd,lnew);
    c = [cupd;cnew];
    a = [aupd;anew];

    % Pruning
%     [r,x,P] = Pruning(r,x,P);
    
    % Recycling
%     [r,x,P,lambdau,xu,Pu] = Recycling(r,x,P,lambdau,xu,Pu);
    
    % Target state extraction
    xest{t} = stateExtract(r_hat,x_hat,model);
    
    % Performance evaluation using GOSPA metric
    [gospa_vals(t,:)] = gospa_dist(get_comps(xlog{t},[1 3]),...
        get_comps(xest{t},[1 3]),gospa_c,gospa_p,gospa_alpha);
    t
    gospa_vals(t,:)
end
