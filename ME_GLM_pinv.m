function [beta_est, sig2_est] = ME_GLM_pinv(Y, X, V, msg)
% _
% Maximum Likelihood Estimation of General Linear Model (using pinv)
% FORMAT [beta_est, sig2_est] = ME_GLM_pinv(Y, X, V, msg)
% 
%     Y   - an n x v data matrix of v time series with n data points
%     X   - an n x p design matrix of p regressors with n data points
%     V   - an n x n covariance matrix embodying covariance assumptions
%     msg - a string used as a message on the SPM progress bar
% 
%     beta_est - a p x v matrix (estimated regression coefficients)
%     sig2_est - a 1 x v vector (estimated residual variance)
% 
% FORMAT [beta_est, sig2_est] = ME_GLM_pinv(Y, X, V) returns the maximum
% likelihood parameter estimates for a general linear model with data Y,
% design matrix X, covariance matrix V.
% 
% The given fMRI data (y) are modelled as a linear combination (b) of
% experimental factors and potential confounds (X), where errors (e) are
% assumed to be normally distributed around zero and to have a known
% covariance structure (V), but unknown residual variance (s2):
%     y = Xb + e, e ~ N(0, s2 V)
% This gives rise to the following likelihood function:
%     p(y|b,s2) = N(y; Xb, s2 V)
% 
% Maximum Likelihood (ML) estimation proceeds by maximizing the (log)
% likelihood function (log) p(y|b,s2) with respect to the free model
% parameters (b,s2).
% 
% When the covariance matrix is equal to the identity matrix (V = I),
% errors are i.i.d. and ML parameter estimates are equivalent to the
% Ordinary Least Squares (OLS) solution [1,2]:
%     beta_est = (X'*X)^-1 * X'*y
%     sig2_est = 1/n * (y-Xb_est)'*(y-Xb_est)
% 
% When the covariance matrix is not equal to the identity matrix (V <> I),
% so that errors are not i.i.d., ML estimates are equivalent to the
% Weighted Least Squares (WLS) solution [1,2]:
%     beta_est = (X'*inv(V)*X)^-1 * X'*inv(V)*y
%     sig2_est = 1/n * (y-Xb_est)'*inv(V)*(y-Xb_est)
% 
% An unbiased estimate for the residual variance is given by:
%     sig2_ub  = n/(n-p) * sig2_est
% 
% References:
% [1] Bishop CM (2006): "Pattern Recognition and Machine Learning".
%     Springer, ch. 3.1, pp. 140-143.
% [2] Koch KR (2007): "Introduction to Bayesian Statistics".
%     Springer, ch. 4.2.2/3, pp. 93-96.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 29/10/2014, 14:05 (V0.2/V7)
%  Last edit: 19/08/2020, 12:04 (V1.4/V20)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 4, msg = 'Estimate GLM ...'; end;
Finter = spm('FigName','ME_GLM_pinv: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P = inv(V);                     % inverse of covariance matrix
W = sqrtm(P);                   % matrix sqrt of precision matrix

% Calculate parameter estimates
%-------------------------------------------------------------------------%
beta_est = pinv(W*X) * (W*Y);   % estimated regression coefficients
resi_est = Y - X*beta_est;      % estimated residuals/errors/noise
sig2_est = zeros(1,v);          % estimated residual variance
for j = 1:v
    sig2_est(j) = 1/n * (W*resi_est(:,j))' * (W*resi_est(:,j));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
% An implementation without a loop would look like this:
%   sig2_est = 1/n * sum(resi_est.^2), if V = I
% However, it would be impossible to make a progress bar then.

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');