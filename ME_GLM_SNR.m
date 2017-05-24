function [mf_SNR, mb_SNR] = ME_GLM_SNR(Y, X, V, B)
% _
% Signal-to-Noise Ratio for Classical General Linear Model
% FORMAT [mf_SNR, mb_SNR] = ME_GLM_SNR(Y, X, V, B)
% 
%     Y - an n x v data matrix of v time series with n data points
%     X - an n x p design matrix of p regressors with n data points
%     V - an n x n covariance matrix embodying covariance assumptions
%     B - a  p x v parameter matrix with WLS regression coefficients
% 
%     mf_SNR - a 1 x v vector of model-free  signal-to-noise ratios
%     mb_SNR - a 1 x v vector of model-based signal-to-noise ratios
% 
% FORMAT [mf_SNR, mb_SNR] = ME_GLM_SNR(Y, X, V, B) calculates model-free
% as well as model-based signal-to-noise ratios for a general linear model
% with data Y, design matrix X and covariance matrix V.
% 
% The model-free signal-to-noise ratio (SNR-mf) is defined as the inverse
% of the coefficient of variation (CV), i.e. as the ratio of mean against
% standard deviation of the measured signal.
% 
% The model-based signal-to-noise ratio (SNR-mb) is defined as the ratio of
% signal variance against noise variance, i.e. as variance of the explained
% signal divided by the variance of the residual signal.
% 
% Please note that the parameter estimates B must be consistent with the
% covariance assumptions V. This means, they have to be calculated as
%     B = (X'*inv(V)*X)^-1 * X'*inv(V)*Y .
% This weighted least squares (WLS) approach is equivalent to whitening
% data and design with a whitening matrix W = sqrtm(inv(V)) before analysis
% and reduces to ordinary least squares (OLS) or maximum likelihood (ML)
% estimation where errors are assumed i.i.d., such that V = eye(n):
%     B = (X'*X)^-1 * X'*Y ,
% If this is case, the input variable V can also be left empty.
% 
% Further information:
%     help ME_GLM
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 11/03/2015, 03:05 (V0.3/V10)
%  Last edit: 11/05/2017, 19:55 (V1.0/V16)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(X,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Whitening of data
%-------------------------------------------------------------------------%
if ~isempty(V)
    P  = inv(V);                % precision matrix
    W  = sqrtm(full(P));        % whitening matrix
    WY = W*Y;                   % whitened data
    WX = W*X;                   % whitened design
else
    P  = eye(n);                % identity matrix
    WY = Y;                     % whitened data
    WX = X;                     % whitened design
end;

% Calculate residuals
%-------------------------------------------------------------------------%
WY_est = WX*B;                 % predicted signal
WE_est = WY - WY_est;          % residual signal

% Calculate model-free signal-to-noise ratio
%-------------------------------------------------------------------------%
mf_SNR = abs(mean(Y))./std(Y);

% Calculate model-based signal-to-noise ratio
%-------------------------------------------------------------------------%
mb_SNR = var(WY_est)./var(WE_est);