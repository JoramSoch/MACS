function [LGBF, post] = ME_BMS_FFX(LME, prior)
% _
% Fixed Effects Bayesian Model Selection using Bayes' Rule
% FORMAT [LGBF, post] = ME_BMS_FFX(LME, prior)
% 
%     LME   - an N x M x S array with log model evidences
%     prior - a  1 x M vector of prior model probabilities
% 
%     LGBF  - an M x M matrix with log group Bayes factors
%     post  - a  1 x M vector of posterior model probabilities
% 
% FORMAT [LGBF, post] = ME_BMS_FFX(LME, prior) returns the posterior model
% probabilities as well as log group Bayes factors, given some prior model
% probabilities and log model evidences for N subjects, M models and S
% sessions.
% 
% Assuming independence across subjects and sessions, log model evidences
% are summed along these dimensions. This gives a likelihood function over
% models which is then used to estimate posterior model probabilities [2].
% 
% The log group Bayes factor of one model against another is the sum of all
% individual log Bayes factors or, in other words, the difference between
% the group log model evidences of these models [1].
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 20/11/2014, 16:50 (V0.2/V8)
%  Last edit: 17/08/2015, 22:00 (V0.3/V11)


% Get model dimensions
%-------------------------------------------------------------------------%
N = size(LME,1);                % number of subjects
M = size(LME,2);                % number of models
S = size(LME,3);                % number of sessions

% Assume independence
%-------------------------------------------------------------------------%
LME = sum(LME,3);               % sum over sessions
LME = sum(LME,1);               % sum over subjects

% Set prior if required
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(prior), prior = 1/M * ones(1,M); end;

% Estimate log group Bayes factors
%-------------------------------------------------------------------------%
LGBF = zeros(M,M);
for j1 = 1:M
    for j2 = 1:M
        LGBF(j1,j2) = LME(j1) - LME(j2);
    end;
end;

% Subtract average log model evidence
%-------------------------------------------------------------------------%
LME = LME - mean(LME);

% Ensure model evidences are computable
%-------------------------------------------------------------------------%
max_val = log(realmax('double'));
for j = 1:M
    LME(j) = sign(LME(j)) * min([max_val,abs(LME(j))]);
end;

% Exponentiate log model evidences
%-------------------------------------------------------------------------%
exp_LME = exp(LME);

% Estimate posterior distribution
%-------------------------------------------------------------------------%
post = exp_LME.*prior;
post = post./sum(post);