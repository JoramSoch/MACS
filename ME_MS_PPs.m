function [post] = ME_MS_PPs(LME, prior)
% _
% Posterior Probabilities within Model Space
% FORMAT [post] = ME_MS_PPs(LME, prior)
% 
%     LME   - an M x V matrix with log model evidences
%     prior - an M x V matrix of prior model probabilities
% 
%     post  - an M x V matrix with posterior probabilities
% 
% FORMAT [post] = ME_MS_PPs(LME, prior) computes posterior model probabilities [1]
% for a number of log model evidences with certain prior model probabilities.
% 
% This procedure simply uses Bayes' theorem to invert probabilities p(y|m),
% i.e. LMEs, into probabilities p(m|y), i.e. PPs. As posterior probabilities
% do not depend on absolute values of LMEs, but only on their relative
% differences, voxel-wise averages are subtracted from LMEs [2].
% 
% Please note that, if the same priors are assumed for voxel, prior can
% also be an M x 1 vector instead of an M x V matrix.
% 
% Further information:
%     help ME_BMS_FFX
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% [2] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/05/2018, 16:10 (V1.2/V18)
%  Last edit: 04/05/2018, 17:10 (V1.2/V18)


% Get family dimensions
%-------------------------------------------------------------------------%
M = size(LME,1);                % number of models
V = size(LME,2);                % number of instances

% Set prior if required
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(prior)
    prior = 1/M * ones(M,1);    % assume uniform prior
end;

% Expand prior if necessary
%-------------------------------------------------------------------------%
if size(prior,2) == 1           % make prior an M x V matrix
    prior = repmat(prior,[1 V]);
end;

% Normalize prior to one
%-------------------------------------------------------------------------%
prior = prior ./ repmat(sum(prior,1),[M 1]);

% Calculate posterior probabilities
%-------------------------------------------------------------------------%
LMEp = LME - repmat(mean(LME,1),[M 1]);
LMEp = exp(LMEp) .* prior;
post = LMEp ./ repmat(sum(LMEp,1),[M 1]);