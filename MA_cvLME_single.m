function MA_cvLME_single(SPM, data, disc, AnC)
% _
% Cross-Validated Log Model Evidence for General Linear Model (single-session)
% FORMAT MA_cvLME_single(SPM, data, disc, AnC)
%     SPM  - a structure specifying an estimated GLM
%     data - a string indicating which data to use (see below)
%     disc - an integer indicating how many volumes to discard (see below)
%     AnC  - a logical indicating accuracy and complexity computation
% 
% FORMAT MA_cvLME_single(SPM, data, disc, AnC) generates a cross-validated
% log model evidence map for a single-session GLM specified by SPM, using
% data indicated by data, discarding a number of volumes indicated by disc
% and calculating accuracy and complexity if AnC is true.
% 
% The present procedure splits the data set into two parts and uses the
% first (second) one to calculate parameter priors for calculating the
% Bayesian log model evidence on the second (first) one. Assumming
% independence between the two parts, the total (cross-validated) log model
% evidence is then equal to the sum of the individual log model evidences.
% 
% The input variable "data" is a string indicating which data to use:
%     If data is 'y',   then the raw data are used.
%     If data is 'Wy',  then the whitened data are used.
%     If data is 'Ky',  then the filtered data are used.
%     If data is 'KWy', then the whitened and filtered data are used.
%     If data is 'WKy', then the filtered and whitened data are used.
% 
% The default for this variable is 'Ky' which means that the analysis
% operates on the filtered data and uses the whitening matrix for
% non-sphericity correction. This is recommended if high-pass filter
% settings are invariant across all models in the model space and
% accomodates for different AR assumptions due to SPM's ReML algorithm.
% 
% The input variable "disc" is an integer indicating how much volumes are
% left out ("discarded") in the middle of the fMRI data set. The default
% for this variable is set such that at least 10 scans are left out and the
% number of remaining scans is divisible by 10, i.e. disc = 10 + mod(n,10).
% 
% The input variable "AnC" is a logical indicating whether model accuracy
% and model complexity are calculated and written to images. The log model
% evidence is the difference of accuracy and complexity: LME = Acc - Com.
% The default for this variable is false.
% 
% Further information:
%     help ME_GLM_NG
%     help ME_GLM_NG_LME
%     help ME_GLM_NG_AnC
% 
% Exemplary usage:
%     MA_cvLME_single(SPM, 'Ky', 15, true);
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469-489.
% [2] Soch J, Meyer AP, Haynes JD, Allefeld C (2017): "How to improve parameter estimates in 
%     GLM-based fMRI data analysis: cross-validated Bayesian model averaging".
%     NeuroImage, vol. 158, pp. 186-195.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/03/2014, 19:00 (V0.1/V1)
%  Last edit: 11/12/2018, 10:15 (V1.3/V19)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    MA_cvLME_single(SPM);
    return
end;

% Estimate model if necessary
%-------------------------------------------------------------------------%
if ~isfield(SPM.xVi,'V')
    SPM_mat = strcat(SPM.swd,'/','SPM.mat');
    MA_GLM_AR_only(SPM_mat); load(SPM_mat);
    MA_cvLME_single(SPM);
    return
end;

% Set data flag if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(data), data = 'Ky'; end;

% Set disc number if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(disc), disc = 10 + mod(size(SPM.xX.X,1),10); end;

% Inactivate AnC if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(AnC), AnC = false; end;

% Call other function if second-level
%-------------------------------------------------------------------------%
if ~isfield(SPM,'Sess')
    disc = mod(size(SPM.xX.X,1),2);
    MA_cvLME_other(SPM,data,disc,AnC);
    return
end;

% Call other function if multi-run
%-------------------------------------------------------------------------%
if numel(SPM.Sess) > 1
    mode = 'N-1';
    MA_cvLME_multi(SPM,data,mode,AnC);
    return
end;

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
Finter = spm('FigName','MA_cvLME_single: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load time series
%-------------------------------------------------------------------------%
Y = MA_load_data(SPM,m_ind);
v = numel(m_ind);


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   P A R T I T I O N                       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (1)');

% Preprocess data if required
%-------------------------------------------------------------------------%
if strcmp(data,'y') | strcmp(data,'Wy')
    if size(K.X0,2) > 1         % Delete the lowest frequency in filter to
        K.X0 = K.X0(:,2:end);   % prevent Ln1/2 from being non-invertible
    end;                        % which can happen due to the partition.
    p = p + size(K.X0,2);       % Design matrix is augmented by filter.
end;                            
if strcmp(data,'y')
  % Y = Y;                      % RAW data are used
    X = [X K.X0];               % design matrix must have filter added
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'Wy')
    Y = W*Y;                    % WHITENED data are used
    X = W*[X K.X0];             % design must have filter and be whitened
    P = speye(n);               % precision is equal to identity matrix
end;
if strcmp(data,'Ky')
    Y = spm_filter(K,Y);        % FILTERED data are used
    X = spm_filter(K,X);        % design matrix must be filtered
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'KWy')
    Y = spm_filter(K,W*Y);      % WHITENED and FILTERED data are used
    X = spm_filter(K,W*X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
end;
if strcmp(data,'WKy')
    Y = W*spm_filter(K,Y);      % FILTERED and WHITENED data are used
    X = W*spm_filter(K,X);      % design must be filtered and whitened
    P = speye(n);               % precision is equal to identity matrix
end;

% Partition data into two parts (1)
%-------------------------------------------------------------------------%
if mod(n-disc,2) == 0           % EVEN number of scans
    s1 = [1:(n-disc)/2];        % leave out d volumes in the middle
    s2 = [((n-disc)/2+disc+1):n];
end;
if mod(n-disc,2) == 1           % UNEVEN number of scans
    s1 = [1:(n-disc-1)/2];      % leave out d volumes and the last one
    s2 = [((n-disc-1)/2+disc+1):(n-1)];
end;

% Partition data into two parts (2)
%-------------------------------------------------------------------------%
if numel(s1) == numel(s2)
    Y1 = Y(s1,:);               % time series
    Y2 = Y(s2,:);
    X1 = X(s1,:);               % design matrix
    X2 = X(s2,:);
    P1 = P(s1,s1);              % precision matrix
    P2 = P(s2,s2);
    n1 = numel(s1);             % data points
    n2 = numel(s2);
end;


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (2)');

% Set (non-informative) priors for both parts
%-------------------------------------------------------------------------%
m0 = zeros(p,1);                % flat Gaussian
L0 = exp(-23)*eye(p);
a0 = 0;                         % Jeffrey's prior
b0 = 0;

% Estimate (informative) posteriors from all data
%-------------------------------------------------------------------------%
[mn, Ln, an, bn] = ME_GLM_NG(Y, X, P, m0, L0, a0, b0, 'Estimate posteriors over both parts 1-2');
clear Y X P

% Estimate posteriors from 1st part (as priors for 2nd part)
%-------------------------------------------------------------------------%
[mn1, Ln1, an1, bn1] = ME_GLM_NG(Y1, X1, P1, m0, L0, a0, b0, 'Estimate 1st part posteriors (as 2nd part priors)');

% Estimate posteriors from 2nd part (as priors for 1st part)
%-------------------------------------------------------------------------%
[mn2, Ln2, an2, bn2] = ME_GLM_NG(Y2, X2, P2, m0, L0, a0, b0, 'Estimate 2nd part posteriors (as 1st part priors)');


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   L O G   M O D E L   E V I D E N C E     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: estimate (3)');

% Preallocate images
%-------------------------------------------------------------------------%
oosLME1 = NaN(size(M));
oosLME2 = NaN(size(M));
if AnC
    oosAcc1 = NaN(size(M));
    oosAcc2 = NaN(size(M));
    oosCom1 = NaN(size(M));
    oosCom2 = NaN(size(M));
end;

% Calculate evidence for 1st part (using estimates from 2nd part)
%-------------------------------------------------------------------------%
oosLME1(m_ind) = ME_GLM_NG_LME(P1, Ln2, an2, bn2, Ln, an, bn);

% Calculate accuracy and complexity for 1st part (using 2nd part)
%-------------------------------------------------------------------------%
if AnC, [oosAcc1(m_ind), oosCom1(m_ind)] = ME_GLM_NG_AnC(X1, P1, mn2, Ln2, an2, bn2, mn, Ln, an, bn, 'Compute accuracy and complexity for 1st part'); end;

% Calculate evidence for 2nd part (using estimates from 1st part)
%-------------------------------------------------------------------------%
oosLME2(m_ind) = ME_GLM_NG_LME(P2, Ln1, an1, bn1, Ln, an, bn);

% Calculate accuracy and complexity for 2nd part (using 1st part)
%-------------------------------------------------------------------------%
if AnC, [oosAcc2(m_ind), oosCom2(m_ind)] = ME_GLM_NG_AnC(X2, P2, mn1, Ln1, an1, bn1, mn, Ln, an, bn, 'Compute accuracy and complexity for 2nd part'); end;

% Calculate total model evidence (assuming independence between parts)
%-------------------------------------------------------------------------%
cvLME = oosLME1 + oosLME2;
if AnC
    cvAcc = oosAcc1 + oosAcc2;
    cvCom = oosCom1 + oosCom2;
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_single: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);

% Write log model evidence
%-------------------------------------------------------------------------%
H.fname   = 'MA_cvLME.nii';
H.descrip = 'MA_cvLME_single: cross-validated log model evidence for general linear model with normal-gamma priors (GLM-NG)';
spm_write_vol(H,reshape(cvLME,m_dim));
SPM.MACS.cvLME = H;
H.fname   = 'MA_cvLME_P1.nii';
H.descrip = 'MA_cvLME_single: log model evidence for 1st part based on priors estimated from 2nd part';
spm_write_vol(H,reshape(oosLME1,m_dim));
SPM.MACS.oosLME(1) = H;
H.fname   = 'MA_cvLME_P2.nii';
H.descrip = 'MA_cvLME_single: log model evidence for 2nd part based on priors estimated from 1st part';
spm_write_vol(H,reshape(oosLME2,m_dim));
SPM.MACS.oosLME(2) = H;

% Write accuracy and complexity
%-------------------------------------------------------------------------%
if AnC
    H.fname   = 'MA_cvAcc.nii';
    H.descrip = 'MA_cvLME_single: cross-validated model accuracy for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvAcc,m_dim));
    SPM.MACS.cvAcc = H;
    H.fname   = 'MA_cvCom.nii';
    H.descrip = 'MA_cvLME_single: cross-validated model complexity for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvCom,m_dim));
    SPM.MACS.cvCom = H;
    H.fname   = 'MA_cvAcc_P1.nii';
    H.descrip = 'MA_cvLME_single: model accuracy for 1st part based on priors estimated from 2nd part';
    spm_write_vol(H,reshape(oosAcc1,m_dim));
    SPM.MACS.oosAcc(1) = H;
    H.fname   = 'MA_cvCom_P1.nii';
    H.descrip = 'MA_cvLME_single: model complexity for 1st part based on priors estimated from 2nd part';
    spm_write_vol(H,reshape(oosCom1,m_dim));
    SPM.MACS.oosCom(1) = H;
    H.fname   = 'MA_cvAcc_P2.nii';
    H.descrip = 'MA_cvLME_single: model accuracy for 2nd part based on priors estimated from 1st part';
    spm_write_vol(H,reshape(oosAcc2,m_dim));
    SPM.MACS.oosAcc(2) = H;
    H.fname   = 'MA_cvCom_P2.nii';
    H.descrip = 'MA_cvLME_single: model complexity for 2nd part based on priors estimated from 1st part';
    spm_write_vol(H,reshape(oosCom2,m_dim));
    SPM.MACS.oosCom(2) = H;
end;

% Save SPM structure
%-------------------------------------------------------------------------%
save(strcat(SPM.swd,'/','SPM.mat'), 'SPM', spm_get_defaults('mat.format'));

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);