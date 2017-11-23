function [r,x,P] = dataAssocNew(wupd,rupd,xupd,Pupd,wnew,rnew,xnew,Pnew)

m = length(wnew)/2;     % m: num of measurements
n = length(wupd)/(m+1); % n: num of pre-existing tracks

H = n*(m+1)+m*2;        % num of single target hypotheses
c = zeros(H,1);         % cost of single target hypotheses

% each new track has two single target hypotheses, one of them has weight
% equal to one, hence zero cost
wupd(wupd==0) = realmin;
c(1:n*(m+1)) = -log(wupd);     % cost of hypotheses of pre-existing tracks
c(n*(m+1)+1:end) = -log(wnew);       % cost of hypotheses of new tracks

% construct binary indicator matrix, with size (n+2m, nm+n+2m)
A = zeros(n+2*m,H);
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

% single target hypotheses updating pre-existing tracks
I = u(1:n*(m+1))==1;
r = rupd(I);
x = xupd(:,I);
P = Pupd(:,:,I);

% single target hypotheses updating new tracks
I = u(n*(m+1)+1:end)==1;
r = [r;rnew(I)];
x = [x xnew(:,I)];
P = cat(3,P,Pnew(:,:,I));

end