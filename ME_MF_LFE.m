function [LFE] = ME_MF_LFE(LME, prior)
% _
% Log Family Evidence for Model Family
% FORMAT [LFE] = ME_MF_LFE(LME, prior)
% 
%     LME   - an M x V matrix with log model evidences
%     prior - a  1 x M vector of prior model probabilities
% 
%     LFE   - a  1 x V vector with log family evidences
% 
% FORMAT [LFE] = ME_MF_LFE(LME, prior) computes the log family evidence for
% a number of log model evidences with certain prior model probabilities.
% 
% Assuming prior and posterior addivity of model probabilities into family
% probabilities as well as a uniform prior within families implies that the
% family evidence is the average of all model evidences belonging to that
% family. However, a sum cannot be logarithmized and log model evidences
% cannot be exponentiated in MATLAB when they are too small.
% 
% This procedure selects the maximum log model evidence in a family as a
% reference point and exponentiates only the difference relative to each
% other log model evidence. If the difference is too large, contribution
% from the respective model will be automatically and rightfully ignored.
% 
% When a non-uniform prior within families is assumed, the computation
% gets more complicated. In this case, the log prior has to be added to
% the log model evidences before the above procedure is performed.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/02/2015, 02:05 (V0.3/V10)
%  Last edit: 17/08/2015, 22:15 (V0.3/V11)


% Get family dimensions
%-------------------------------------------------------------------------%
M = size(LME,1);                % number of models
V = size(LME,2);                % number of instances

% Set prior if required
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(prior)
    prior = 1/M * ones(1,M);
end;

% Normalize prior to one
%-------------------------------------------------------------------------%
if sum(prior) ~= 1
    prior = 1/sum(prior) * prior;
end;

% If prior is uniform
%-------------------------------------------------------------------------%
if sum(abs(diff(prior))) == 0
    LME_fam  = LME;
    LME_max  = max(LME_fam,[],1);
    LME_diff = LME_fam - repmat(LME_max,[M 1]);
    LFE      = LME_max + log(mean(exp(LME_diff),1));
end;

% If prior is non-uniform
%-------------------------------------------------------------------------%
if sum(abs(diff(prior))) > 0
    LME = LME + repmat(log(prior)',[1 V]) + log(M);
    LFE = ME_MF_LFE(LME, (1/M)*ones(1,M));
end;