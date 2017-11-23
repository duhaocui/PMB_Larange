clc;clear
dbstop if error

% Set simulation parameters
Pd = 0.5; % probability of detection
lfai = 10; % expected number of false alarms per scan
numtruth = 2; % number of targets
simcasenum = 1; % simulation case 1 or 2 (see paper)
%simcasenum = 2;
if (simcasenum == 1) % covariance used for mid-point initialisation
    Pmid = 1e-6*eye(4);
else
    Pmid = 0.25*eye(4);
end

% Generate truth
[model,measlog,xlog] = gentruth(Pd,lfai,numtruth,Pmid,simcasenum);

% Initialise filter parameters
stateDimensions = size(model.xb,1);
% Multi-Bernoulli representation
n = 0;
r = zeros(n,1);
x = zeros(stateDimensions,n);
P = zeros(stateDimensions,stateDimensions,n);
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
    % Predict
    [r,x,P,lambdau,xu,Pu] = predictStep(r,x,P,lambdau,xu,Pu,model);
    
    % Predict all single target hypotheses of previous scan
    
    % Update
    [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew] = ...
        updateStepNew(lambdau,xu,Pu,r,x,P,measlog{t},model);
    
    % Update all predicted single target hypotheses of previous scan
    
    % Store single target hypotheses of current scan, send to next scan
    % Order: missed detection, single target hypotheses updated by the
    % first, second, third,..., measurement, followed by single target
    % hypothese of newly detected targets

    % Data assciation
    [r,x,P] = dataAssocNew(wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew);

    % Target state extraction
    xest{t} = stateExtract(r,x,model);
    
    % Pruning
    [r,x,P] = Pruning(r,x,P);
    
    % Recycling
    [r,x,P,lambdau,xu,Pu] = Recycling(r,x,P,lambdau,xu,Pu);
    
    % Performance evaluation using GOSPA metric
    [gospa_vals(t,:)] = gospa_dist(get_comps(xlog{t},[1 3]),...
        get_comps(xest{t},[1 3]),gospa_c,gospa_p,gospa_alpha);
end
