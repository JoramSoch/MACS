function [beta_est, sig2_est] = ME_GLM(Y, X, V, msg)
% _
% Estimation of Classical General Linear Model
% FORMAT [beta_est, sig2_est] = ME_GLM(Y, X, V, msg)
% 
%     Y   - an n x v data matrix of v time series with n data points
%     X   - an n x p design matrix of p regressors with n data points
%     V   - an n x n covariance matrix embodying covariance assumptions
%     msg - a string used as a message on the SPM progress bar
% 
%     beta_est - a p x v matrix (estimated regression coefficients)
%     sig2_est - a 1 x v vector (estimated residual variance)
% 
% FORMAT [beta_est, sig2_est] = ME_GLM(Y, X, V) returns the classical
% parameter estimates for a general linear model with data Y, design
% matrix X, covariance matrix V.
% 
% The given fMRI data (y) are modelled as a linear combination (b) of
% experimental factors and potential confounds (X), where errors (e) are
% assumed to be normally distributed around zero and to have a known
% covariance structure (V):
%     y = Xb + e, e ~ N(0, V)
% This gives rise to the following likelihood function:
%     p(y|beta) = N(y; Xb, V)
% This implies the following residual sum of squares:
%     RSS(beta) = (y-Xb)'(y-Xb) = |y-Xb|^2
% 
% Maximum Likelihood (ML) or Ordinary Least Squares (OLS) estimation
% proceeds by minimizing RSS with respect to beta which is equivalent
% to the residual vector being orthogonal to the design space (X'e = 0).
% ML/OLS estimates are given by:
%     beta_est = (X'*X)^-1 * X'*y
%     sig2_est = 1/n * (y-Xb_est)'*(y-Xb_est)
% 
% When the covariance matrix is not equal to the identity matrix (V <> I),
% so that errors are not i.i.d., Gauss-Markov (GM) estimation or Weighted
% Least Squares (WLS) must be used. GM/WLS estimates are given by:
%     beta_est = (X'*inv(V)*X)^-1 * X'*inv(V)*y
%     sig2_est = 1/n * (y-Xb_est)'*inv(V)*(y-Xb_est)
% 
% An unbiased estimate for the residual variance is given by:
%     sig2_est = 1/(n-p) * (y-Xb_est)'*(y-Xb_est)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 29/10/2014, 14:05 (V0.2/V7)
%  Last edit: 10/03/2015, 19:00 (V0.3/V10)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 4, msg = 'Estimate GLM ...'; end;
Finter = spm('FigName','ME_GLM: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P = inv(V);                     % inverse of data covariance matrix

% Estimate posterior parameters
%-------------------------------------------------------------------------%
beta_est = (X'*P*X)^-1 * X'*P*Y;% estimated regression coefficients
resi_est = Y - X*beta_est;      % estimated residuals/errors/noise
sig2_est = zeros(1,v);          % estimated residual variance
for j = 1:v
    sig2_est(j) = 1/n * (resi_est(:,j))' * P * (resi_est(:,j));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
% An implementation without a loop would look like this:
%   sig2_est = 1/n * sum(resi_est.^2), if V = I
% However, it would be impossible to make a progress bar then.

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');