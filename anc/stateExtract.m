function [est] = stateExtract(r,x,model)
% Extract target states

ss = r > model.existThresh;
est = x(:,ss);

end

