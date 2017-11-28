function [est] = stateExtract(r,x,model)
% Extract target states

% ss = r > model.existThresh;
% est = x(:,ss);

% MAP cardinality estimate
r(r==1) = 1-eps;                    % remove numerical error
ss = false(size(r));
pcard = prod(1-r)*poly(-r./(1-r));
[~,n] = max(pcard);
[~,o] = sort(-r);
n = n - 1;
ss(o(1:n)) = true;
est = x(:,ss);


end

