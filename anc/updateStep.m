function [lambdau,xu,Pu,wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew] = updateStep(lambdau,xu,Pu,r,x,P,z,model)
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
%Output:
% lambdau(k), xu(:,k) and Pu(:,:,k) give the updated intensity, state
%  estimate and covariance for the k-th mixture component of the unknown
%  target Poisson Point Process (PPP)

% Extract parameters from model
Pd = model.Pd;
H = model.H;
R = model.R;
lambda_fa = model.lambda_fa;
lambdab_threshold = 1e-3;

% Interpret sizes from inputs
n = length(r);
[stateDimensions,nu] = size(xu);
[measDimensions,m] = size(z);

% Allocate memory for existing tracks
% Note no gating--assume any measurement can go with any track
% Obviously this would be done differently in practice
wupd = zeros(n,m+1);
rupd = zeros(n,m+1);
xupd = zeros(stateDimensions,n,m+1);
Pupd = zeros(stateDimensions,stateDimensions,n,m+1);

% Allocate memory for new tracks
wnew = zeros(m,1);
rnew = zeros(m,1);
xnew = zeros(stateDimensions,m);
Pnew = zeros(stateDimensions,stateDimensions,m);

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
    wupd(i,1) = 1 - r(i) + r(i)*(1-Pd);
    rupd(i,1) = r(i)*(1-Pd)/wupd(i,1);
    xupd(:,i,1) = x(:,i);
    Pupd(:,:,i,1) = P(:,:,i);
    
    % Create hypotheses with measurement updates
    S = H*P(:,:,i)*H' + R;
    sqrt_det2piS = sqrt(det(2*pi*S));
    K = P(:,:,i)*H'/S;
    Pplus = P(:,:,i) - K*H*P(:,:,i);
    for j = 1:m
        v = z(:,j) - H*x(:,i);
        wupd(i,j+1) = r(i)*Pd*exp(-0.5*v'/S*v)/sqrt_det2piS;
        rupd(i,j+1) = 1;
        xupd(:,i,j+1) = x(:,i) + K*v;
        Pupd(:,:,i,j+1) = Pplus;
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
    %i = n + j;
    for k = 1:nu
        v = z(:,j) - H*xu(:,k);
        ck(k) = lambdau(k)*Pd*exp(-0.5*v'/Sk(:,:,k)*v)/sqrt_det2piSk(k);
        yk(:,k) = xu(:,k) + Kk(:,:,k)*v;
    end
    C = sum(ck);
    wnew(j) = C + lambda_fa;
    rnew(j) = C/wnew(j);
    ck = ck/C;
    xnew(:,j) = yk*ck;
    for k = 1:nu
        v = xnew(:,j) - yk(:,k);
        Pnew(:,:,j) = Pnew(:,:,j) + ck(k)*(Pk(:,:,k) + v*v');
    end
end

% Update (i.e., thin) intensity of unknown targets
lambdau = (1-Pd)*lambdau;

% Not shown in paper--truncate low weight components
ss = lambdau > lambdab_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);
