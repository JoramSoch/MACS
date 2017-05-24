function [PLL, MAPLL] = ME_GLM_DIC(Y, X, V)
% _
% Deviance Information Criterion for the General Linear Model
% FORMAT [PLL, MAPLL] = ME_GLM_DIC(Y, X, V)
% 
%     Y  - an n x v data matrix of v time series with n data points
%     X  - an n x p design matrix of p regressors with n data points
%     V  - an n x n covariance matrix embodying temporal correlations
% 
%     PLL   - a 1 x v vector of posterior log-likelihood values
%     MAPLL - a 1 x v vector of maximum-a-posteriori log-likelihoods
% 
% FORMAT [PLL, MAPLL] = ME_GLM_DIC(Y, X, V) calculates components of the
% deviance information criterion (DIC) for the general linear model (GLM)
% with data Y, design matrix X and covariance matrix V.
% 
% The DIC is given by
%     DIC = -2*MAPLL + 2*pD
% where pD is given by
%     pD  = -2*PLL + 2*MAPLL
% where PLL is the posterior log-likelihood, i.e. the posterior expec-
% tation of the log-likelihood, and MAPLL is the maximum-a-posteriori
% log-likelihood, i.e. the log-likelihood at the posterior modes.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 11/05/2017, 23:20 (V1.0/V16)
%  Last edit: 11/05/2017, 23:20 (V1.0/V16)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters
P = inv(V);

% Set non-informative prior
%-------------------------------------------------------------------------%
m0 = zeros(p,1);                % flat Gaussian
L0 = zeros(p,p);
a0 = 0;                         % Jeffrey's prior
b0 = 0;

% Estimate Bayesian GLM
%-------------------------------------------------------------------------%
[mn, Ln, an, bn] = ME_GLM_NG(Y, X, P, m0, L0, a0, b0, 'Estimate Bayesian GLM for DIC calculation');

% Calculate posterior log-likelihood
%-------------------------------------------------------------------------%
[PLL, Com] = ME_GLM_NG_AnC(X, P, m0, L0, a0, b0, mn, Ln, an, bn, 'Calculate posterior log-likelihood');
clear Com

% Calculate maximum-a-posteriori log-likelihood
%-------------------------------------------------------------------------%
beta_MAP = mn;
tau_MAP  = (an-1)./bn;
MAPLL    = ME_GLM_LL(Y, X, V, beta_MAP, 1./tau_MAP);
clear beta_MAP tau_MAP