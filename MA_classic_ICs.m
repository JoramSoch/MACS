function MA_classic_ICs(SPM, data, ICs)
% _
% Classical Information Criteria for General Linear Model
% FORMAT MA_classic_ICs(SPM, data, ICs)
%     SPM  - a structure specifying an estimated GLM
%     data - a string indicating which data to use (see below)
%     ICs  - a cell array indicating information criteria to calculate
% 
% FORMAT MA_classic_ICs(SPM, data, ICs) calculates classical information
% criteria that only depend on the maximum log-likelihood (MLL), the number
% of data points (n) and the number of model parameters (p).
% 
% The input variable "data" is a string indicating which data to use:
%     If data is 'y',   then the raw data are used.
%     If data is 'Wy',  then the whitened data are used.
%     If data is 'Ky',  then the filtered data are used.
%     If data is 'KWy', then the whitened and filtered data are used.
%     If data is 'WKy', then the filtered and whitened data are used.
% 
% The input variable "ICs" is a cell array of strings and can contain:
%     If IC is 'AIC',  then it's Akaike information criterion [1].
%     If IC is 'AICc', then it's corrected Akaike information criterion [2].
%     If IC is 'BIC',  then it's Bayesian information criterion [3].
%     If IC is 'DIC',  then it's Deviance information criterion [4].
%     If IC is 'HQC',  then it's Hannan-Quinn information criterion [5].
%     If IC is 'KIC',  then it's Kullback information criterion [6].
%     If IC is 'KICc', then it's corrected Kullback information criterion [7].
% 
% Note that the DIC is - in contrast to the BIC - a "real" Bayesian
% information criterion that requires prior distributions. Although it
% therefore depends on more parameters than MLL, n and p, it is included
% in this routine for convenience and because it can be written in a
% similar form as AIC and BIC.
% 
% Further information:
%     help ME_GLM_LL
%     help ME_GLM_GoF
%     help ME_GLM_SNR
% 
% Exemplary usage:
%     MA_classic_ICs(SPM, 'Ky', {'AIC', 'BIC'});
% 
% References:
% [1] Akaike H (1974): "A new look at the statistical model identification".
%     IEEE Transactions on Automatic Control, vol. 19, pp. 716-723.
% [2] Hurvich CM, Tsai CL (1989): "Regression and Time Series Model Selection
%     in Small Samples". Biometrika, vol. 76, pp. 297-307.
% [3] Schwarz G (1978): "Estimating the Dimension of a Model".
%     The Annals of Statistics, vol. 6, pp. 461-464.
% [4] Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A (2002): 
%     "Bayesian measures of model complexity and fit".
%     Journal of the Royal Statistical Society B, vol. 64, pp. 583-639.
% [5] Hannan EJ, Quinn BG (1979):
%     "The Determination of the order of an autoregression".
%     Journal of the Royal Statistical Society B, vol. 41, pp. 190-195.
% [6] Cavanaugh JE (1999): "A large-sample model selection
%     criterion based on Kullback's symmetric divergence".
%     Statistics and Probability Letters, vol. 42, p. 333-344.
% [7] Cavanaugh JE (2004): "Criteria for linear model selection
%     based on Kullback's symmetric divergence".
%     Australian & New Zealand Journal of Statistics, vol. 46, pp. 257-274.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 09:30 (V0.99/V15)
%  Last edit: 09/03/2018, 14:05 (V1.2/V18)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    MA_classic_ICs(SPM);
    return
end;

% Estimate model if necessary
%-------------------------------------------------------------------------%
if ~isfield(SPM.xVi,'V')
    SPM = spm_spm(SPM);
    MA_classic_ICs(SPM);
    return
end;

% Set data flag if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(data), data = 'Ky'; end;

% Set criterion if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(ICs), ICs = {'AIC', 'BIC'}; end;

% Change to SPM.swd if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Get model parameters
%-------------------------------------------------------------------------%
X = SPM.xX.X;                   % design matrix
K = SPM.xX.K;                   % filtering matrix
W = SPM.xX.W;                   % whitening matrix
V = SPM.xVi.V;                  % non-sphericity
n = size(X,1);                  % number of observations
p = size(X,2);                  % number of regressors

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_classic_ICs: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load time series
%-------------------------------------------------------------------------%
Y = MA_load_data(SPM,m_ind);
v = length(m_ind);

% Load parameter estimates
%-------------------------------------------------------------------------%
B = MA_load_betas(SPM,m_ind);


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_classic_ICs: estimate');

% Preprocess data if required
%-------------------------------------------------------------------------%
if strcmp(data,'y')
  % Y = Y;                      % RAW data are used
    X = [X K.X0];               % design matrix must have filter added
    P = spm_inv(V);             % precision is inverse of non-sphericity
  % V = V;                      % covariance remains at non-sphericity
end;
if strcmp(data,'Wy')
    Y = W*Y;                    % WHITENED data are used
    X = W*[X K.X0];             % design must have filter and be whitened
    P = speye(n);               % precision is equal to identity matrix
    V = speye(n);               % covariance is equal to identity matrix
end;
if strcmp(data,'Ky')
    Y = spm_filter(K,Y);        % FILTERED data are used
    X = spm_filter(K,X);        % design matrix must be filtered
    P = spm_inv(V);             % precision is inverse of non-sphericity
  % V = V;                      % covariance remains at non-sphericity
end;
if strcmp(data,'KWy')
    Y = spm_filter(K,W*Y);      % WHITENED and FILTERED data are used
    X = spm_filter(K,W*X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
    V = speye(n);               % covariance is equal to identity matrix    
end;
if strcmp(data,'WKy')
    Y = W*spm_filter(K,Y);      % FILTERED and WHITENED data are used
    X = W*spm_filter(K,X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
    V = speye(n);               % covariance is equal to identity matrix
end;

% Calculate maximum log-likelihood
%-------------------------------------------------------------------------%
if ~isempty(find(strcmp(ICs','DIC')==0))
    [s2, R2, adj_R2, gen_R2] = ME_GLM_GoF(Y, X, V, B);
    MLL = ME_GLM_LL(Y, X, V, B, s2);
    p   = p + 1;                % The number of model parameters also
    clear R2 adj_R2 gen_R2      % includes the residual variance!
end;

% Calculate posterior log-likelihood
%-------------------------------------------------------------------------%
if ~isempty(cell2mat(strfind(ICs','DIC')))
    [PLL, LLP] = ME_GLM_DIC(Y, X, V);
    pD = -2*PLL + 2*LLP;
    clear PLL
end;

% Calculate information criteria
%-------------------------------------------------------------------------%
map = cell(numel(ICs),1);
for i = 1:numel(ICs)
    map{i} = NaN(size(M));
    switch ICs{i}
        case 'AIC'              % Akaike information criterion
            XIC = -2*MLL + 2*p;
        case 'AICc'             % corrected Akaike information criterion
            XIC = -2*MLL + 2*p + (2*p.*(p+1))./(n-p-1);
        case 'BIC'              % Bayesian information criterion
            XIC = -2*MLL + p.*log(n);
        case 'DIC'              % Deviance information criterion
            XIC = -2*LLP + 2*pD;
        case 'HQC'              % Hannan-Quinn information criterion
            XIC = -2*MLL + 2*p.*log(log(n));
        case 'KIC'              % Kullback information criterion
            XIC = -2*MLL + 3*p;
        case 'KICc'             % corrected Kullback information criterion
            XIC = -2*MLL + n.*log(n./(n-p)) + (n.*((n-p).*(2*p+3)-2))/((n-p-2).*(n-p));
    end;
    map{i}(m_ind) = XIC;
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_classic_ICs: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);

% Write information criterion
%-------------------------------------------------------------------------%
for i = 1:numel(ICs)
    H.fname   = strcat('MA_ICs_',ICs{i},'.nii');
    H.descrip = sprintf('MA_classic_ICs: %s',ICs{i});
    spm_write_vol(H,reshape(map{i},m_dim));
    eval(strcat('SPM.MACS.',ICs{i},' = H;'));
end;

% Save SPM structure
%-------------------------------------------------------------------------%
save(strcat(SPM.swd,'/','SPM.mat'), 'SPM', spm_get_defaults('mat.format'));

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);