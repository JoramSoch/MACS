function [mn, Ln, an, bn] = ME_GLM_NG(Y, X, P, m0, L0, a0, b0, msg)
% _
% Estimation of General Linear Model with Normal-Gamma Priors
% FORMAT [mn, Ln, an, bn] = ME_GLM_NG(Y, X, P, m0, L0, a0, b0, msg)
% 
%     Y   - an n x v data matrix of v time series with n data points
%     X   - an n x p design matrix of p regressors with n data points
%     P   - an n x n precision matrix embodying covariance assumptions
%     m0  - a  p x v matrix (prior mean for regression coefficients)
%     L0  - a  p x p matrix (prior precision for regression coefficients)
%     a0  - a  1 x 1 scalar (prior shape for residual variance)
%     b0  - a  1 x v vector (prior rate for residual variance)
%     msg - a string used as a message on the SPM progress bar
% 
%     mn  - a  p x v matrix (posterior mean for regression coefficients)
%     Ln  - a  p x p matrix (posterior precision for regression coefficients)
%     an  - a  1 x 1 scalar (posterior shape for residual variance)
%     bn  - a  1 x v vector (posterior rate for residual variance)
% 
% FORMAT [mn, Ln, an, bn] = ME_GLM_NG(Y, X, P, m0, L0, a0, b0) returns the
% posterior parameter estimates for a general linear model with data Y,
% design matrix X, precision matrix P and normal-gamma distributed priors
% for regression coefficients (m0, L0) and residual variance (a0, b0).
% 
% Please note that, if the same priors are assumed for each time series,
% m0 can also be a p x 1 vector and b0 can also be a 1 x 1 scalar.
% 
% The given fMRI data (y) are modelled as a linear combination (b) of
% experimental factors and potential confounds (X), where errors (e) are
% assumed to be normally distributed around zero and to have a known
% covariance structure (V), but unknown variance factor (sigma^2):
%     y = Xb + e, e ~ N(0, sigma^2*V)
% For mathematical convenience, we write:
%     sigma^2 = 1/tau, V = P^–1 => sigma^2*V = (tau*P)^–1
% This gives rise to the following likelihood function:
%     p(y|beta,tau) = N(y; Xb, (tau*P)^–1)
% We use normal-gamma distributions as conjugate priors [1,2]:
%     p(beta|tau) = N(beta; m0, (tau*L0)^–1)
%     p(tau) = Gam(tau; a0, b0)
% This results in normal-gamma distributions for the posteriors:
%     p(beta|tau,y) = N(beta; mn, (tau*Ln)^–1)
%     p(tau|y) = Gam(tau; an, bn)
% Parameter posteriors as well as log model evidence are analytic [3].
% 
% References:
% [1] Bishop CM (2006): "Pattern Recognition and Machine Learning".
%     Springer, ch. 3.3, pp. 175-177.
% [2] Koch KR (2000): "Einführung in die Bayes-Statistik".
%     Springer, ch. 4.3.2, pp. 117-119.
% [3] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 472–473, eq. 6/9.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 16:30 (V0.2/V6)
%  Last edit: 15/09/2016, 07:10 (V0.9a/V13a)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Enlarge priors if required
%-------------------------------------------------------------------------%
if size(m0,2) == 1              % make m0 a p x v matrix
    m0 = repmat(m0,[1 v]);
end;
if size(b0,2) == 1              % make b0 a 1 x v vector
    b0 = b0*ones(1,v);
end;

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 8, msg = 'Estimate GLM ...'; end;
Finter = spm('FigName','ME_GLM_NG: estimate');
spm_progress_bar('Init', 100, msg, '');

% Estimate posterior parameters
%-------------------------------------------------------------------------%
Ln = X'*P*X + L0;               % precision for regression coefficients
mn = inv(Ln) * (X'*P*Y + L0*m0);% mean for regression coefficients
an = a0 + n/2;                  % shape for residual variance
bn = zeros(1,v);                % rate for residual variance
for j = 1:v
    bn(j) = b0(j) + 1/2*(Y(:,j)'*P*Y(:,j) + m0(:,j)'*L0*m0(:,j) - mn(:,j)'*Ln*mn(:,j));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
% An implementation without a loop would look like this:
%   bn = b0 + 1/2* (sum(Y.*(P*Y)) + sum(m0.*(L0*m0)) - sum(mn.*(Ln*mn)));
% However, it would be impossible to make a progress bar then.

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');