function [sig2, R2, adj_R2, gen_R2] = ME_GLM_GoF(Y, X, V, B)
% _
% Goodness of Fit for Classical General Linear Model
% FORMAT [sig2, R2, adj_R2, gen_R2] = ME_GLM_GoF(Y, X, V, B)
% 
%     Y - an n x v data matrix of v time series with n data points
%     X - an n x p design matrix of p regressors with n data points
%     V - an n x n covariance matrix embodying covariance assumptions
%     B - a  p x v parameter matrix with WLS regression coefficients
% 
%       sig2 - a 1 x v vector of residual variances
%         R2 - a 1 x v vector of coefficients of determination
%     adj_R2 - a 1 x v vector of adjusted R^2 values
%     gen_R2 - a 1 x v vector of generalized R^2 values
% 
% FORMAT [sig2, R2, adj_R2, gen_R2] = ME_GLM_GoF(Y, X, V, B) calculates
% residual variance sig2, the coefficient of determination R^2 as well as
% adjusted R^2 and generalized R^2 for a general linear model with data Y,
% design matrix X and covariance matrix V.
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
% First edit: 04/03/2015, 15:15 (V0.3/V10)
%  Last edit: 18/08/2017, 16:50 (V1.1/V17)


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
WY_est  = WX*B;                 % predicted signal
WE_est  = WY - WY_est;          % residual signal
WY_mean = mean(WY);             % average signal

% Init progress bar
%-------------------------------------------------------------------------%
msg = 'Compute sums of squares ...';
Finter = spm('FigName','ME_GLM_GoF: estimate');
spm_progress_bar('Init', 100, msg, '');

% Calculate sums of squares
%-------------------------------------------------------------------------%
sig2   = zeros(1,v);            % residual variance
SS_tot = zeros(1,v);            % total sum of squares
SS_res = zeros(1,v);            % residual sum of squares
SS_reg = zeros(1,v);            % regressed sum of squares
for j = 1:v
    sig2(j)   = 1/n * sum(WE_est(:,j).^2);
    SS_res(j) =   n * sig2(j); 
  % SS_res(j) = sum((WE_est(:,j)).^2);
    SS_tot(j) = sum((WY(:,j) - WY_mean(j)).^2);
  % SS_reg(j) = sum((WY_est(:,j) - WY_mean(j)).^2);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
% An implementation without a loop would look like this:
%   sig2   = 1/n * sum(WE_est.^2);
%   SS_res = sum((WE_est).^2); % = n * sig2
%   SS_tot = sum((WY - repmat(WY_mean,[n 1])).^2);
%   SS_reg = sum((WY_est - repmat(WY_mean,[n 1])).^2);
% However, it would be impossible to make a progress bar then.

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Calculate R^2 and adj. R^2
%-------------------------------------------------------------------------%
    R2 = 1 - (SS_res)./(SS_tot); % = SS_reg./SS_tot
adj_R2 = 1 - (SS_res./(n-p)) ./ (SS_tot./(n-1));

% Init progress bar
%-------------------------------------------------------------------------%
% msg = 'Compute log-likelihoods ...';
% Finter = spm('FigName','ME_GLM_GoF: estimate');
% spm_progress_bar('Init', 100, msg, '');

% Calculate log-likelihoods
%-------------------------------------------------------------------------%
% LL0 = zeros(1,v);
% LLm = zeros(1,v);
% for j = 1:v
%     LL0(j) = - n/2 * log(var(Y(:,j),1)) - 1/2 * (1/var(Y(:,j),1)) * ...
%               (Y(:,j) - mean(Y(:,j)))' * P * (Y(:,j) - mean(Y(:,j)));
%     LLm(j) = - n/2 * log(sig2(j)) - 1/2 * (1/sig2(j)) * ...
%               (Y(:,j) - X*B(:,j))' * P * (Y(:,j) - X*B(:,j));
%     if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
% end;
% LL0 = LL0 - 1/2 * log(det(V)) - n/2 * log(2*pi);
% LLm = LLm - 1/2 * log(det(V)) - n/2 * log(2*pi);

% Clear progress bar
%-------------------------------------------------------------------------%
% spm_progress_bar('Clear');

% Calculate generalized R^2
%-------------------------------------------------------------------------%
gen_R2 = R2; % = 1 - exp(LL0 - LLm).^(2/n)
% In a classical GLM, generalized R^2 is identical with R^2.