function [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew] = updateStepNew(lambdau,xu,Pu,r,x,P,z,model)
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

% Allocate memory for new tracks
wnew = zeros(2*m,1);
rnew = zeros(2*m,1);
xnew = zeros(stateDimensions,2*m);
Pnew = zeros(stateDimensions,stateDimensions,2*m);

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
end

% Update (i.e., thin) intensity of unknown targets
lambdau = (1-Pd)*lambdau;

% Not shown in paper--truncate low weight components
ss = lambdau > lambdab_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);
