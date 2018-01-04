function [lambdau,xu,Pu,rupd,xupd,Pupd,lupd,cupd,rnew,xnew,Pnew,lnew,cnew,aupd,anew] = ...
    updateStepNew2(lambdau,xu,Pu,r,x,P,l,c,z,a,model)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENTS
%Syntax: [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew] =
%          update(lambdau,xu,Pu,r,x,P,z,model)
%Input:
% lambdau(k), xu(:,k) and Pu(:,:,k) give the intensity, state estimate and
%  covariance for the k-th mixture component of the unknown target Poisson
%  Point Process (PPP)
% r(i), x(:,i) and P(:,:,i) give the probability of existence, state
%   estimate and covariance for the i-th multi-Bernoulli component (track)
% z(:,j) is measurement j
% model is a structure containing parameters of measurement model

% Extract parameters from model
Pd = model.Pd;
H = model.H;
R = model.R;
lambda_fa = model.lambda_fa;
lambdab_threshold = 1e-4;

% Interpret sizes from inputs
% If it is a hypothesis with zero existence probability, no update
% e.g. n_invalid, n_valid, nupd = n_invalid + n_valid*(m+1);
[stateDimensions,nu] = size(xu);
[measDimensions,m] = size(z);
n = length(r);  % number of single target hypotheses to be updated
% length(r)=length(find(r~=0))+length(find(r==0))
nupd = length(find(r~=0))*m + n; 

% Allocate memory for existing tracks (single target hypotheses)
% Note no gating--assume any measurement can go with any track
% Obviously this would be done differently in practice
wupd = zeros(nupd,1);
rupd = zeros(nupd,1);
xupd = zeros(stateDimensions,nupd);
Pupd = zeros(stateDimensions,stateDimensions,nupd);
lupd = cell(nupd,1);
cupd = zeros(nupd,1);
% Keep record ancestor information, i.e., the index of the single target 
% hypothesis being updated, can be used in N-scan pruning
% Currently, only the latest ancestor information is recorded, i.e., 1-scan
aupd = zeros(nupd,1);

% Implement algorithm

% Update existing tracks
iupd = 0;   % initiate updating index
for i = 1:n
    % first determine whether the single target hypothesis has valid
    % existence probability or not
    if r(i) == 0
        iupd = iupd+1;
        wupd(iupd) = 1;
        cupd(iupd) = 0;
        rupd(iupd) = 0;
        xupd(:,iupd) = zeros(stateDimensions,1);
        Pupd(:,:,iupd) = zeros(stateDimensions,stateDimensions,1);
        lupd{iupd} = [l{i};0];
        aupd(iupd) = a(i);
    else
        % Create missed detection hypothesis
        iupd = iupd+1;
        wupd(iupd) = 1 - r(i) + r(i)*(1-Pd);
        cupd(iupd) = c(i)-log(wupd(iupd));
        rupd(iupd) = r(i)*(1-Pd)/wupd(iupd);
        xupd(:,iupd) = x(:,i);
        Pupd(:,:,iupd) = P(:,:,i);
        % If it is missed detection, add label 0
        lupd{iupd} = [l{i};0];
        aupd(iupd) = a(i);
        
        % Create hypotheses with measurement updates
        S = H*P(:,:,i)*H' + R;
        sqrt_det2piS = sqrt(det(2*pi*S));
        K = P(:,:,i)*H'/S;
        Pplus = P(:,:,i) - K*H*P(:,:,i);
        for j = 1:m
            iupd = iupd+1;
            v = z(:,j) - H*x(:,i);
            wupd(iupd) = r(i)*Pd*exp(-0.5*v'/S*v)/sqrt_det2piS;
            cupd(iupd) = c(i)-log(wupd(iupd));
            rupd(iupd) = 1;
            xupd(:,iupd) = x(:,i) + K*v;
            Pupd(:,:,iupd) = Pplus;
            % Otherwise, add the index of the measurement at current scan
            lupd{iupd} = [l{i};j];
            aupd(iupd) = a(i);
        end
    end
end

% Prune single target hypothesis with really small likelihood
idx_keep = wupd>1e-6;
rupd = rupd(idx_keep);
xupd = xupd(:,idx_keep);
Pupd = Pupd(:,:,idx_keep);
lupd = lupd(idx_keep);
cupd = cupd(idx_keep);
aupd = aupd(idx_keep);

% Allocate memory for new tracks, each new track contains two single target
% hypothese
wnew = zeros(2*m,1);
rnew = zeros(2*m,1);
xnew = zeros(stateDimensions,2*m);
Pnew = zeros(stateDimensions,stateDimensions,2*m);
lnew = cell(2*m,1);
anew = zeros(2*m,1);

% Allocate temporary working for new tracks
Sk = zeros(measDimensions,measDimensions,nu);
Kk = zeros(stateDimensions,measDimensions,nu);
Pk = zeros(stateDimensions,stateDimensions,nu);
ck = zeros(nu,1);
sqrt_det2piSk = zeros(nu,1);
yk = zeros(stateDimensions,nu);

% Create a new track for each measurement by updating PPP with measurement
for k = 1:nu
    Sk(:,:,k) = H*Pu(:,:,k)*H' + R;
    sqrt_det2piSk(k) = sqrt(det(2*pi*Sk(:,:,k)));
    Kk(:,:,k) = Pu(:,:,k)*H'/Sk(:,:,k);
    Pk(:,:,k) = Pu(:,:,k) - Kk(:,:,k)*H*Pu(:,:,k);
end
for j = 1:m
    for k = 1:nu
        v = z(:,j) - H*xu(:,k);
        ck(k) = lambdau(k)*Pd*exp(-0.5*v'/Sk(:,:,k)*v)/sqrt_det2piSk(k);
        yk(:,k) = xu(:,k) + Kk(:,:,k)*v;
    end
    C = sum(ck);
    % first single target hypothesis for measurement associated to previous
    % track, second for new track
    wnew(2*j-1) = 1;
    wnew(2*j) = C + lambda_fa;
    rnew(2*j) = C/wnew(2*j);
    ck = ck/C;
    xnew(:,2*j) = yk*ck;
    for k = 1:nu
        v = xnew(:,2*j) - yk(:,k);
        Pnew(:,:,2*j) = Pnew(:,:,2*j) + ck(k)*(Pk(:,:,k) + v*v');
    end
    % for trajectory purpose
    lnew{2*j-1} = 0;    % add 0, if there is no new target
    lnew{2*j} = j;      % otherwise, add measurement index
    if isempty(aupd)
        anew(2*j-1:2*j) = j;
    else  
        anew(2*j-1:2*j) = aupd(end)+j;
    end
end
cnew = -log(wnew);

% If there is no pre-existing track, no need to create non-exist hypothesis
% if isempty(rupd)
%     rnew = rnew(2:2:end);
%     xnew = xnew(:,2:2:end);
%     Pnew = Pnew(:,:,2:2:end);
%     lnew = lnew(2:2:end);
%     cnew = cnew(2:2:end);
%     anew = anew(2:2:end);
% end

% Update (i.e., thin) intensity of unknown targets
lambdau = (1-Pd)*lambdau;

% Truncate low weight components
ss = lambdau > lambdab_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);
