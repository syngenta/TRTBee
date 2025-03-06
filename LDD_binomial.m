function [LD50, slope] = LDD_binomial(x, yInit, yEnd)
% Function to fit log-likelihood using binomial method
%% Inputs:
% x
% yInit: initial survivors
% yEnd: survivors at the end of the experiment
%%
p = @(b,x) 1./(1 + exp(b(1)*(log(x) - log(b(2))))); % LDD equation

logLik = @(b) -sum( nchoosek_binomial_fitting(yInit, yEnd) + yEnd.*log(p(b,x)) + (yInit - yEnd).*log(1 - p(b,x)));
% yInit could be 1 if normalized
% p is the test value from the log-log function
% logLik = @(b) sum(log(nchoosek(1, y)) + yEnd.*log(p) + (1 - y).*log(1 - p));

% Options for optimisation algorithm
opts = optimset('MaxFunEvals',500000, 'MaxIter',1000000);
startVals = [-1, min(x)];

B = fminsearch(logLik, startVals, opts); % run algorithm


LD50 = B(2); slope = B(1); % Extract results
