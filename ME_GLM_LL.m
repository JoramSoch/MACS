function LL = ME_GLM_LL(Y, X, V, B, s2)
% _
% Log-Likelihood Function for the General Linear Model
% FORMAT LL = ME_GLM_LL(Y, X, V, B, s2)
% 
%     Y  - an n x v data matrix of v time series with n data points
%     X  - an n x p design matrix of p regressors with n data points
%     V  - an n x n covariance matrix embodying temporal correlations
%     B  - a  p x v matrix of regression coefficients
%     s2 - a  1 x v vector of residual variances
% 
%     LL - a  1 x v vector of log-likelihood values
% 
% FORMAT LL = ME_GLM_LL(Y, X, V, B, s2) calculates the log-likelihood
% for the general linear model with data Y, design matrix X, covariance
% matrix V, regression coefficients B and residual variances s2.
% 
% Further information:
%     help ME_GLM
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
% First edit: 18/03/2017, 08:40 (V0.99/V15)
%  Last edit: 26/04/2019, 17:50 (V1.4/V20)


% Number of data sets
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
d = floor(v/100);

% Precision and residuals
%-------------------------------------------------------------------------%
P = inv(V);                     % precision = inverse of covariance
R = Y - X*B;                    % residuals = measurement minus prediction

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ME_GLM_LL: estimate');
spm_progress_bar('Init', 100, 'Calculate log-likelihood...', '');

% Log-likelihood function
%-------------------------------------------------------------------------%
LL_con = 1/2 * MD_mvn_logdet(P,true) - n/2 * log(2*pi*s2);
LL_var = zeros(1,v);
for j = 1:v
    LL_var(j) = - 1/2 * ( R(:,j)' * P * R(:,j) ) / s2(j);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
LL = LL_con + LL_var;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');