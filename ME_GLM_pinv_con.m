function [h, p, stats] = ME_GLM_pinv_con(Y, X, V, c, type, alpha, msg)
% _
% Contrast-Based Inference for Classical General Linear Model (using pinv)
% FORMAT [h, p, stats] = ME_GLM_pinv_con(Y, X, V, c, type, alpha, msg)
% 
%     h     - a 1 x v vector of test results (H0: false; H1: true)
%     p     - a 1 x v vector of p-values for the t- or F-test
%     stats - a structure variable with the following fields:
%     o tstat - value of the t-statistic  (if type is 't')
%     o df    - t-test degrees of freedom (if type is 't')
%     o Fstat - value of the F-statistic       (if type is 'F')
%     o df_n  - numerator degrees of freedom   (if type is 'F')
%     o df_d  - denominator degrees of freedom (if type is 'F')
% 
%     Y     - an n x v data matrix (n data points, v variables)
%     X     - an n x p design matrix (n data points, p regressors)
%     V     - an n x n correlation matrix (n x n observations)
%     c     - a  p x 1 contrast vector or a p x q contrast matrix
%     type  - a  string indicating the type of inference ('t' or 'F')
%     alpha - a  scalar indicating the significance level (e.g. 0.05)
%     msg   - a  string used as a message on the SPM progress bar
% 
% FORMAT [h, p, stats] = ME_GLM_pinv_con(Y, X, V, c, type, alpha, msg)
% performns a statistical test indicated by type and specified by the
% contrast c using the design matrix X, the estimated model parameters
% beta_est and sig2_est and returns test result h, p-value p and stats,
% information about the test statistic.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 29/10/2018, 15:55 (V1.3/V19)
%  Last edit: 19/08/2020, 12:14 (V1.4/V20)


% Modify contrast if necessary
%-------------------------------------------------------------------------%
if numel(c) > 1 && size(c,1) == 1
    c = c';                     % transform to column vector
end;
if size(c,1) ~= size(X,2)       % append zeros to contrast
    c = [c; zeros(size(X,2)-size(c,1),size(c,2))];
end;

% Specify type if necessary
%-------------------------------------------------------------------------%
if nargin < 5 || isempty(type)
    if size(c,2) == 1           % contrast vector
        type = 't';             % t-test
    else                        % contrast matrix
        type = 'F';             % F-test
    end;
end;

% Specify alpha if necessary
%-------------------------------------------------------------------------%
if nargin < 6 || isempty(alpha)
    alpha = 0.05;               % significance level
end;

% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(X,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Estimate model parameters
%-------------------------------------------------------------------------%
[beta_est, sig2_est] = ME_GLM_pinv(Y, X, V);
 sig2_ub = (n/(n-p)) * sig2_est;

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 7, msg = 'Calculate contrasts ...'; end;
Finter = spm('FigName','ME_GLM_pinv_con: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare inference
%-------------------------------------------------------------------------%
P    = inv(V);                  % precision matrix (inverse of covariance)
W    = sqrtm(P);                % whitening matrix (matrix sqrt of precision)
covB = pinv(W*X) * pinv(W*X)';  % covariance matrix (of the beta estimates)

% Perform t-test
%-------------------------------------------------------------------------%
if strcmp(type,'t')
    c_cov_c = c'*covB*c;        % quadratic form
    con_est = c'*beta_est;      % contrast estimate
    den_est = zeros(1,v);       % denominator value
    for j = 1:v
        den_est(j) = sqrt(sig2_ub(j) * c_cov_c);
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
        % The loop is not necessary, but I want a progress bar.
    end;
    stats.tstat = con_est./den_est;
    stats.df    = n-p;
    p = 1 - spm_Tcdf(stats.tstat, stats.df);
    h = p < alpha;
end;

% Perform F-test
%-------------------------------------------------------------------------%
if strcmp(type,'F')
    C = c;
    q = size(C,2);
    C_cov_C = C'*covB*C;        % quadratic form
    inv_CcC = inv(C_cov_C);     % inverse of that
    con_est = C'*beta_est;      % contrast estimate
    num_est = zeros(1,v);       % numerator value
    for j = 1:v
        num_est(j) = con_est(:,j)' * inv_CcC * con_est(:,j);
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
        % The loop is not necessary, but I want a progress bar.
    end;
    stats.Fstat = (1/q) * num_est./sig2_ub;
    stats.df_n  = q;
    stats.df_d  = n-p;
    p = 1 - spm_Fcdf(stats.Fstat, stats.df_n, stats.df_d);
    h = p < alpha;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');