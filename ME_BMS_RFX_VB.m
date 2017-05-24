function [alpha_post] = ME_BMS_RFX_VB(LME, alpha_prior)
% _
% Random Effects Bayesian Model Selection using Variational Bayes
% FORMAT [alpha_post] = ME_BMS_RFX_VB(LME, alpha_prior)
% 
%     LME         - an N x M x S array with log model evidences
%     alpha_prior - a  1 x M vector of prior Dirichlet parameters
% 
%     alpha_post  - a  1 x M vector of posterior Dirichlet parameters
% 
% FORMAT [alpha_post] = ME_BMS_RFX_VB(LME, alpha_prior) estimates a random-
% effects (RFX) Bayesian model selection (BMS) using Variational Bayes (VB)
% and returns the Dirichlet parameters of the posterior distribution over
% models. These parameters describe how often the model has been observed
% in the population.
% 
% The structure of the model is given as follows [1,2,3]:
%     p(y|m) = 1st level model evidences as 2nd level likelihood function
%     p(m|r) = Mult(m; 1, r) where r = [r1, ..., rM] model probabilities
%     p(r|a) = Dir(r; alpha) where a = [a1, ..., aM] Dirichlet parameters
%     alpha0 = [1, ..., 1] giving a uniform prior distribution over models
% 
% Model estimation proceeds using Variational Bayes. Essentially, this
% script implements the procedure described in [1,2] and should give the
% same results as spm_BMS.m while being substantially faster [4].
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Rosa MJ, Bestmann S, Harrison L, Penny W (2010):
%     "Bayesian model selection maps for group studies".
%     NeuroImage, vol. 49, pp. 217-224.
% [3] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [4] Soch J & Allefeld C (2015): "Non-Critical Comments on
%     Bayesian Model Selection". Internal Report, June 2015.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 20/11/2014, 17:10 (V0.2/V8)
%  Last edit: 15/09/2016, 09:05 (V0.9a/V13a)


% Get model dimensions
%-------------------------------------------------------------------------%
N = size(LME,1);                % number of subjects
M = size(LME,2);                % number of models
S = size(LME,3);                % number of sessions

% Assume independence
%-------------------------------------------------------------------------%
LME = sum(LME,3);               % sum over sessions

% Set prior if required
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(alpha_prior), alpha_prior = ones(1,M); end;

% Subtract average log model evidence
%-------------------------------------------------------------------------%
LME = LME - repmat(mean(LME,2),[1 M]);

% Ensure model evidences are computable
%-------------------------------------------------------------------------%
max_val = log(realmax('double'));
for i = 1:N
    for j = 1:M
        LME(i,j) = sign(LME(i,j)) * min([max_val,abs(LME(i,j))]);
    end;
end;

% Initialize model estimation
%-------------------------------------------------------------------------%
alpha = alpha_prior;
conv  = 1e-4;
diff  = 1;

% Variational Bayes algorithm
%-------------------------------------------------------------------------%
while diff > conv
    % compute u_nm
    u_nm = exp( LME + repmat(psi(alpha),[N 1]) - psi(sum(alpha))*ones(N,M) );
    % compute g_nm
    g_nm = u_nm ./ repmat(sum(u_nm,2),[1 M]);
    % update alpha
    alpha_old = alpha;
    alpha = alpha_prior + sum(g_nm,1);
    diff = norm(alpha - alpha_old);
end;

% Get posterior distribution
%-------------------------------------------------------------------------%
alpha_post = alpha;