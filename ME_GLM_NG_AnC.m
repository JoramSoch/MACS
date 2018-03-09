function [Acc, Com] = ME_GLM_NG_AnC(X, P, m0, L0, a0, b0, mn, Ln, an, bn, msg)
% _
% Accuracy and Complexity for General Linear Model with Normal-Gamma Priors
% FORMAT [Acc, Com] = ME_GLM_NG_AnC(X, P, m0, L0, a0, b0, mn, Ln, an, bn, msg)
% 
%     X   - an n x p design matrix of p regressors with n data points
%     P   - an n x n precision matrix embodying covariance assumptions
%     m0  - a  p x v matrix (prior mean for regression coefficients)
%     L0  - a  p x p matrix (prior precision for regression coefficients)
%     a0  - a  1 x 1 scalar (prior shape for residual variance)
%     b0  - a  1 x v vector (prior rate for residual variance)
%     mn  - a  p x v matrix (posterior mean for regression coefficients)
%     Ln  - a  p x p matrix (posterior precision for regression coefficients)
%     an  - a  1 x 1 scalar (posterior shape for residual variance)
%     bn  - a  1 x v vector (posterior rate for residual variance)
%     msg - a string used as a message on the SPM progress bar
% 
%     Acc - a  1 x v vector of (Bayesian) model accuracies
%     Com - a  1 x v vector of (Bayesian) model complexities
% 
% FORMAT [Acc, Com] = ME_GLM_NG_AnC(X, P, m0, L0, a0, b0, mn, Ln, an, bn, msg)
% returns model accuracy and model complexity [1] for a general linear model
% with design matrix X, precision matrix P and normal-gamma distributed
% priors/posteriors for regression coefficients (m0, L0 / mn, Ln) and
% residual variance (a0, b0 / an, bn).
% 
% Further information:
%     help ME_GLM_NG
%     help ME_GLM_NG_LME
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469-489, eqs. C.2/C.4.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/01/2015, 14:10 (V0.3/V9)
%  Last edit: 09/03/2018, 12:15 (V1.2/V18)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(mn,2);                 % number of time series
n = size(X,1);                  % number of data points
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
if nargin < 11, msg = 'Compute accuracy and complexity ...'; end;
Finter = spm('FigName','ME_GLM_NG_AnC: estimate');
spm_progress_bar('Init', 100, msg, '');

% Calculate parameter error
%-------------------------------------------------------------------------%
PE = zeros(1,v);
for j = 1:v
    PE(j) = (m0(:,j)-mn(:,j))' * L0 * (m0(:,j)-mn(:,j));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;

% Calculate accuracy
%-------------------------------------------------------------------------%
Acc = - 1/2 * (an./bn) .* (2*(bn-b0) - PE) + ...
        1/2 * log(det(P)) - 1/2 * trace(X'*P*X * inv(Ln)) - ...
        n/2 * log(2*pi) + n/2 * (psi(an)-log(bn));

% Calculate complexity
%-------------------------------------------------------------------------%
Com = + 1/2 * (an./bn) .* (PE - 2*(bn-b0)) - ...
        1/2 * log(det(L0)/det(Ln)) + 1/2 * trace(L0*inv(Ln)) - p/2 + ...
         a0 * log(bn./b0) - (gammaln(an)-gammaln(a0)) + (an-a0)*psi(an);

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');