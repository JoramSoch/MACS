function [alpha_post, exp_freq, exc_prob] = ME_BMS_RFX_GS(LME, alpha_prior, N_samp)
% _
% Random Effects Bayesian Model Selection using Gibbs Sampling
% FORMAT [alpha_post, exp_freq, exc_prob] = ME_BMS_RFX_GS(LME, alpha_prior, N_samp)
% 
%     LME         - an N x M x S array with log model evidences
%     alpha_prior - a  1 x M vector of prior Dirichlet parameters
%     N_samp      - an integer indicating the number of samples
% 
%     alpha_post  - a  1 x M vector of posterior Dirichlet parameters
%     exp_prob    - a  1 x M vector of expected frequencies over models
%     exc_prob    - a  1 x M vector of exceedance probabilities over models
% 
% FORMAT [alpha_post, exp_freq, exc_prob] = ME_BMS_RFX_GS(LME, alpha_prior, N_samp)
% performs random-effects (RFX) Bayesian model selection (BMS) by Gibbs Sampling (GS)
% and returns the Dirichlet parameters of the posterior distribution over
% models as well as expected frequencies and exceedance probabilities over
% models. These probabilities describe how strong the data favor a model in
% comparison with all the others.
% 
% The structure of the model is given as follows [1,2,3]:
%     p(y|m) = 1st level model evidences as 2nd level likelihood function
%     p(m|r) = Mult(m; 1, r) where r = [r1, ..., rM] model probabilities
%     p(r|a) = Dir(r; alpha) where a = [a1, ..., aM] Dirichlet parameters
%     alpha0 = [1, ..., 1] giving a uniform prior distribution over models
% 
% Model estimation proceeds using Gibbs Sampling. Essentially, this script
% implements the procedure described in [3] and should give the same
% results as spm_BMS_gibbs.m.
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
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 20/11/2014, 19:10 (V0.2/V8)
%  Last edit: 17/08/2015, 22:00 (V0.3/V11)


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

% Set number of samples if required
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(N_samp), N_samp = 2e4; end;

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
alpha = zeros(N_samp,M);
r     = zeros(N_samp+1,M);
r(1,:)= MD_dirrnd(alpha_prior,1);

% Gibbs Sampling algorithm
%-------------------------------------------------------------------------%
for k = 1:N_samp
    
    % Sample m given r and y
    %---------------------------------------------------------------------%
    % compute u_nm
    u_nm = exp( LME + repmat(log(r(k,:)),[N 1]) );
    % compute g_nm
    g_nm = u_nm ./ repmat(sum(u_nm,2),[1 M]);
    % sample m|r,y
    a_nm = zeros(N,M);
    for i = 1:N
        j = MD_multrnd(g_nm(i,:),1);
        a_nm(i,j) = 1;
    end;
    
    % Sample r given m and y
    %---------------------------------------------------------------------%
    % compute beta
    beta = sum(a_nm,1);
    % update alpha
    alpha(k,:) = alpha_prior + beta;
    % sample r|m,y
    r(k+1,:) = MD_dirrnd(alpha(k,:),1);
    
end;

% Get posterior distribution
%-------------------------------------------------------------------------%
alpha_post = mean(alpha(floor(N_samp/2):end,:));

% Estimate expected probability
%-------------------------------------------------------------------------%
exp_freq = mean(r(floor(N_samp/2):end,:));

% Estimate exceedance probability
%-------------------------------------------------------------------------%
r_samp   = r(floor(N_samp/2):end,:);
[mr, jr] = max(r_samp,[],2);
exc_prob = zeros(1,M);
for j = 1:M
    exc_prob(j) = sum(jr==j)./size(r_samp,1);
end;