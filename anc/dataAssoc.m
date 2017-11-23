function [r,x,P] = dataAssoc(wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew)

[n,m] = size(wupd);     % n: num of pre-existing tracks
m = m-1;                % m: num of measurements

H = n*(m+1)+m*2;        % num of single target hypotheses
c = zeros(H,1);         % cost of single target hypotheses

% each new track has two single target hypotheses, one of them has weight
% equal to one, hence zero cost
c(1:n*(m+1)) = -log(reshape(wupd'+realmin,[n*(m+1),1]));     % cost of hypotheses of pre-existing tracks
c(n*(m+1)+2:2:end) = -log(wnew);                             % cost of hypotheses of new tracks

% construct binary indicator matrix, with size (n+2m, nm+n+2m)
A = zeros(n+2*m,n*(m+1)+2*m);
% for each pre-existing track
for i = 1:n
    A(i,(i-1)*(m+1)+1:i*(m+1)) = 1;
end
% for each new track
for i = 1:m
    A(n+i,n*(m+1)+(i-1)*2+1:n*(m+1)+i*2) = 1;
end
% for each measurement
for i = 1:m
    A(n+m+i,i+1:m+1:n*(m+1)) = 1;       % for each pre-existing track
    A(n+m+i,n*(m+1)+2*i) = 1;           % for each new track
end

% constriant equality
b = ones(n+2*m,1);

% binary indicator vector, decides which single target hypotheses are
% included in the "best" global hypotheses. 
options = optimoptions('intlinprog','Display','off');
u = intlinprog(c,1:length(c),[],[],A,b,zeros(length(c),1),ones(length(c),1),[],options);

% selected single target hypotheses
stateDimention = size(xnew,1);
r = zeros(0,1);
x = zeros(stateDimention,0);
P = zeros(stateDimention,stateDimention,0);

% single target hypotheses updating pre-existing tracks
Apre = u(1:n*(m+1));
for i = 1:n
    I = Apre((i-1)*(m+1)+1:i*(m+1))==1;
    r = [r;rupd(i,I)];
    x = [x xupd(:,i,I)];
    P = cat(3,P,Pupd(:,:,i,I));
end

% single target hypotheses updating new tracks
Anew = u(n*(m+1)+1:end);
for i = 1:m
    I = Anew(2*i-1:2*i)==1;
    if I(2) == true
        r = [r;rnew(i)];
        x = [x xnew(:,i)];
        P = cat(3,P,Pnew(:,:,i));
    end
end

end