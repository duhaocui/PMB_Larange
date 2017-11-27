function [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,lupd,aupd,cupd,wnew,rnew,xnew,Pnew,lnew,cnew] = ...
    updateStepNew(lambdau,xu,Pu,r,x,P,l,c,z,model)
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
n = length(r);
[stateDimensions,nu] = size(xu);
[measDimensions,m] = size(z);

% Allocate memory for existing tracks
% Note no gating--assume any measurement can go with any track
% Obviously this would be done differently in practice
wupd = zeros(n*(m+1),1);
rupd = zeros(n*(m+1),1);
xupd = zeros(stateDimensions,n*(m+1));
Pupd = zeros(stateDimensions,stateDimensions,n*(m+1));
lupd = cell(n*(m+1),1);
cupd = zeros(n*(m+1),1);
% Keep record ancestor information, used in N-scan pruning
% Currently, only the latest ancestor information is recorded
aupd = zeros(n*(m+1),1);

% Allocate memory for new tracks
wnew = zeros(2*m,1);
rnew = zeros(2*m,1);
xnew = zeros(stateDimensions,2*m);
Pnew = zeros(stateDimensions,stateDimensions,2*m);
lnew = cell(2*m,1);
cnew = zeros(2*m,1);

% No need to record new track ancestor
% anew = zeros(2*m,1);

% Allocate temporary working for new tracks
Sk = zeros(measDimensions,measDimensions,nu);
Kk = zeros(stateDimensions,measDimensions,nu);
Pk = zeros(stateDimensions,stateDimensions,nu);
ck = zeros(nu,1);
sqrt_det2piSk = zeros(nu,1);
yk = zeros(stateDimensions,nu);

% Implement algorithm

% Update existing tracks
for i = 1:n
    % Create missed detection hypothesis
    wupd((i-1)*m+i) = 1 - r(i) + r(i)*(1-Pd);
    rupd((i-1)*m+i) = r(i)*(1-Pd)/wupd((i-1)*m+i);
    xupd(:,(i-1)*m+i) = x(:,i);
    Pupd(:,:,(i-1)*m+i) = P(:,:,i);
    % If it is missed detection, add label 0
    lupd{(i-1)*m+i} = [l{i};0];
    aupd((i-1)*m+i) = i;
    cupd((i-1)*m+i) = c(i)-log(wupd((i-1)*m+i));
    
    % Create hypotheses with measurement updates
    S = H*P(:,:,i)*H' + R;
    sqrt_det2piS = sqrt(det(2*pi*S));
    K = P(:,:,i)*H'/S;
    Pplus = P(:,:,i) - K*H*P(:,:,i);
    for j = 1:m
        v = z(:,j) - H*x(:,i);
        wupd(m*i-m+i+j) = r(i)*Pd*exp(-0.5*v'/S*v)/sqrt_det2piS;
        rupd(m*i-m+i+j) = 1;
        xupd(:,m*i-m+i+j) = x(:,i) + K*v;
        Pupd(:,:,m*i-m+i+j) = Pplus;
        % Otherwise, add the index of the measurement at current scan
        lupd{m*i-m+i+j} = [l{i};j];
        aupd(m*i-m+i+j) = i;
        
        if wupd(m*i-m+i+j)==0
            wupd(m*i-m+i+j) = realmin;
        end
        cupd(m*i-m+i+j) = c(i)-log(wupd(m*i-m+i+j));
    end
    
end

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
end
cnew = -log(wnew);

% Update (i.e., thin) intensity of unknown targets
lambdau = (1-Pd)*lambdau;

% Not shown in paper--truncate low weight components
ss = lambdau > lambdab_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);
