function MA_cvLME_multi(SPM, data, mode, AnC)
% _
% Cross-Validated Log Model Evidence for General Linear Model (multi-session)
% FORMAT MA_cvLME_multi(SPM, data, mode, AnC)
%     SPM  - a structure specifying an estimated GLM
%     data - a string indicating which data to use (see below)
%     mode - a string indicating how to cross-validate (see below)
%     AnC  - a logical indicating accuracy and complexity computation
% 
% FORMAT MA_cvLME_multi(SPM, data, mode, AnC) generates a cross-validated
% log model evidence map for a multi-session GLM specified by SPM, using
% data indicated by data and using a cross-validation mode indicated by
% mode.
% 
% The present procedure performs N-fold cross-validation and uses N-1
% sessions (or 1 session) to calculate parameter priors that are used for
% calculating the Bayesian log model evidence in the N-th session (or in
% the next session). Assumming independence between sessions, the total
% (cross-validated) log model evidence is then equal to the sum of the
% individual log model evidences.
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
% The input variable "mode" is a string indicating how to cross-validate:
%     If mode is 'N-1',  then N-1 sessions are used to calculate the prior.
%     If mode is 'prev', then the previous session is used to calculate it.
% 
% The default for this variable is 'N-1' which means that all but one
% session are used to calculate parameter posteriors. These posteriors are
% then used as priors to calculate a log model evidence in the remaining
% session. This is repeated for all sessions and log model evidences are
% summed across sessions.
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
%     MA_cvLME_multi(SPM, 'Ky', 'N-1', false);
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% [2] Soch J, Meyer AP, Haynes JD, Allefeld C (2017): "How to improve parameter estimates in 
%     GLM-based fMRI data analysis: cross-validated Bayesian model averaging".
%     NeuroImage, in review. URL: http://biorxiv.org/content/early/2016/12/20/095778
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/03/2014, 19:00 (V0.1/V1)
%  Last edit: 19/02/2018, 12:15 (V1.2/V18)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    MA_cvLME_multi(SPM);
    return
end;

% Estimate model if necessary
%-------------------------------------------------------------------------%
if ~isfield(SPM.xVi,'V')
    SPM_mat = strcat(SPM.swd,'/','SPM.mat');
    MA_GLM_AR_only(SPM_mat); load(SPM_mat);
    MA_cvLME_multi(SPM);
    return
end;

% Set data flag if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(data), data = 'Ky'; end;

% Set CV mode if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(mode), mode = 'N-1'; end;

% Inactivate AnC if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(AnC), AnC = false; end;

% Call other function if GLM for EEG
%-------------------------------------------------------------------------%
if ~isfield(SPM,'Sess')
    SPM.Sess = 1;
    disc = 10 + mod(size(SPM.xX.X,1),10);
    MA_cvLME_single(SPM,data,disc,AnC);
    return
end;

% Call other function if single-run
%-------------------------------------------------------------------------%
if length(SPM.Sess) == 1
    disc = 10 + mod(size(SPM.xX.X,1),10);
    MA_cvLME_single(SPM,data,disc,AnC);
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
s = length(SPM.Sess);           % number of sessions
n = zeros(1,s);                 % numbers of volumes
p = zeros(1,s);                 % numbers of parameters
S = struct([]);                 % session structure

% Allocate sessions
%-------------------------------------------------------------------------%
for h = 1:s
    n(h)   = numel(SPM.Sess(h).row);
    p(h)   = numel(SPM.Sess(h).col)+1;
    S(h).n = [SPM.Sess(h).row];
    S(h).p = [SPM.Sess(h).col SPM.xX.iB(h)];
end;

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_multi: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load time series
%-------------------------------------------------------------------------%
Y = MA_load_data(SPM,m_ind);
v = length(m_ind);


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   P A R T I T I O N                       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_multi: estimate (1)');

% Preprocess data if required
%-------------------------------------------------------------------------%
if strcmp(data,'y') | strcmp(data,'Wy')
    X_DCT = blkdiag(K.X0);      % Form discrete cosine transform (DCT) set
    for h = 1:s                 % Add DCT regressors to list of parameters
        p_DCT = size(K(h).X0,2);
        p(h)  = p(h) + p_DCT;
        if h == 1, S(h).p = [S(h).p (size(X,2)+[1:p_DCT])]; end;
        if h  > 1, S(h).p = [S(h).p (max(S(h-1).p)+[1:p_DCT])]; end;
    end;
end;                            
if strcmp(data,'y')
  % Y = Y;                      % RAW data are used
    X = [X X_DCT];              % design matrix must have filter added
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'Wy')
    Y = W*Y;                    % WHITENED data are used
    X = W*[X X_DCT];            % design must have filter and be whitened
    P = speye(sum(n));          % precision is equal to identity matrix
end;
if strcmp(data,'Ky')
    Y = spm_filter(K,Y);        % FILTERED data are used
    X = spm_filter(K,X);        % design matrix must be filtered
    P = spm_inv(V);             % precision is inverse of non-sphericity
end;
if strcmp(data,'KWy')
    Y = spm_filter(K,W*Y);      % WHITENED and FILTERED data are used
    X = spm_filter(K,W*X);      % design must be whitened and filtered
    P = speye(sum(n));          % precision is equal to identity matrix
end;
if strcmp(data,'WKy')
    Y = W*spm_filter(K,Y);      % FILTERED and WHITENED data are used
    X = W*spm_filter(K,X);      % design must be filtered and whitened
    P = speye(sum(n));          % precision is equal to identity matrix
end;

% Partition data into sessions
%-------------------------------------------------------------------------%
for h = 1:s
    S(h).Y = Y(S(h).n,:);       % time series
    S(h).X = X(S(h).n,S(h).p);  % design matrix
    S(h).P = P(S(h).n,S(h).n);  % precision matrix
end;


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_multi: estimate (2)');

% Set non-informative priors
%-------------------------------------------------------------------------%
m0_ni = zeros(p(1),1);          % flat Gaussian
L0_ni = exp(-23)*eye(p(1));
a0_ni = 0;                      % Jeffrey's prior
b0_ni = 0;

% Estimate posteriors from all data
%-------------------------------------------------------------------------%
if strcmp(mode,'N-1')
    X = vertcat(S(1:s).X);
    [mn_all, Ln_all, an_all, bn_all] = ME_GLM_NG(Y, X, P, m0_ni, L0_ni, a0_ni, b0_ni, sprintf('Estimate posteriors over all sessions 1-%d',s));
end;
clear Y X P

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for g = 1:s
    
    % List sessions for this fold
    %---------------------------------------------------------------------%
    fold = 1:s;
    fold = fold(fold~=g);

    % Concatenate design matrices
    %---------------------------------------------------------------------%
    Yf = vertcat(S(fold).Y);
    Xf = vertcat(S(fold).X);
    Pf = blkdiag(S(fold).P);
    
    % If cross-validation mode is "N-1 sessions"
    %---------------------------------------------------------------------%
    if strcmp(mode,'N-1')

        % Set non-informative priors
        %-----------------------------------------------------------------%
        m0 = m0_ni; L0 = L0_ni; a0 = a0_ni; b0 = b0_ni;
        
        % Estimate posteriors for this fold
        %-----------------------------------------------------------------%
        [mn, Ln, an, bn] = ME_GLM_NG(Yf, Xf, Pf, m0, L0, a0, b0, sprintf('Estimate posteriors that serve as priors for session %d',g));

        % Save priors and posteriors for this session
        %-----------------------------------------------------------------%
        S(g).m0 = mn;     S(g).L0 = Ln;     S(g).a0 = an;     S(g).b0 = bn;
        S(g).mn = mn_all; S(g).Ln = Ln_all; S(g).an = an_all; S(g).bn = bn_all;
        clear m0 L0 a0 b0 mn Ln an bn

    end;
    
    % If cross-validation mode is "previous session"
    %---------------------------------------------------------------------%
    if strcmp(mode,'prev')
        
        % Set non-informative priors
        %-----------------------------------------------------------------%
        m0 = m0_ni; L0 = L0_ni; a0 = a0_ni; b0 = b0_ni;

        % Estimate posteriors for previous session
        %-----------------------------------------------------------------%
        h = g - 1; if g == 1, h = s; end;
        [mn, Ln, an, bn] = ME_GLM_NG(S(h).Y, S(h).X, S(h).P, m0, L0, a0, b0, sprintf('Estimate priors for session %d as posteriors from session %d',g,h));

        % Set these as informative priors
        %-----------------------------------------------------------------%
        m0 = mn; L0 = Ln; a0 = an; b0 = bn;

        % Estimate posteriors for current session
        %-----------------------------------------------------------------%
        [mn, Ln, an, bn] = ME_GLM_NG(S(g).Y, S(g).X, S(g).P, m0, L0, a0, b0, sprintf('Estimate posteriors for session %d using priors from session %d',g,h));
        
        % Save priors and posteriors for this session
        %-----------------------------------------------------------------%
        S(g).m0 = m0; S(g).L0 = L0; S(g).a0 = a0; S(g).b0 = b0;
        S(g).mn = mn; S(g).Ln = Ln; S(g).an = an; S(g).bn = bn;
        clear m0 L0 a0 b0 mn Ln an bn
        
    end;
    
end;


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   L O G   M O D E L   E V I D E N C E     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_multi: estimate (3)');

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for g = 1:s
    
    % Preallocate log model evidence for this session
    %---------------------------------------------------------------------%
    oosLME = NaN(size(M));
    if AnC
        oosAcc = NaN(size(M));
        oosCom = NaN(size(M));
    end;
    
    % Calculate log model evidence for this session
    %---------------------------------------------------------------------%
    oosLME(m_ind) = ME_GLM_NG_LME(S(g).P, S(g).L0, S(g).a0, S(g).b0, S(g).Ln, S(g).an, S(g).bn);
    if AnC
        [oosAcc(m_ind), oosCom(m_ind)] = ME_GLM_NG_AnC(S(g).X, S(g).P, S(g).m0, S(g).L0, S(g).a0, S(g).b0, S(g).mn, S(g).Ln, S(g).an, S(g).bn, sprintf('Compute accuracy and complexity for session %d',g));
    end;
    
    % Save log model evidence for this session
    %---------------------------------------------------------------------%
    S(g).LME = oosLME;
    if AnC
        S(g).Acc = oosAcc;
        S(g).Com = oosCom;
    end;
    clear oosLME oosAcc oosCom
    
    % Sum log model evidences from all sessions
    %---------------------------------------------------------------------%
    if g == 1
        cvLME = 0;
        cvAcc = 0;
        cvCom = 0;
    end;
        cvLME = cvLME + S(g).LME;
    if AnC
        cvAcc = cvAcc + S(g).Acc;
        cvCom = cvCom + S(g).Com;
    end;

end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_cvLME_multi: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);

% Write log model evidence
%-------------------------------------------------------------------------%
H.fname   = 'MA_cvLME.nii';
H.descrip = 'MA_cvLME_multi: cross-validated log model evidence for general linear model with normal-gamma priors (GLM-NG)';
spm_write_vol(H,reshape(cvLME,m_dim));
SPM.MACS.cvLME = H;
for g = 1:s
    H.fname   = strcat('MA_cvLME_','S',int2str(g),'.nii');
    H.descrip = sprintf('MA_cvLME_multi: log model evidence for session %s based on priors estimated from other sessions',int2str(g));
    spm_write_vol(H,reshape(S(g).LME,m_dim));
    SPM.MACS.oosLME(g) = H;
end;

% Write accuracy and complexity
%-------------------------------------------------------------------------%
if AnC
    H.fname   = 'MA_cvAcc.nii';
    H.descrip = 'MA_cvLME_multi: cross-validated model accuracy for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvAcc,m_dim));
    SPM.MACS.cvAcc = H;
    H.fname   = 'MA_cvCom.nii';
    H.descrip = 'MA_cvLME_multi: cross-validated model complexity for general linear model with normal-gamma priors (GLM-NG)';
    spm_write_vol(H,reshape(cvCom,m_dim));
    SPM.MACS.cvCom = H;
    for g = 1:s
        H.fname   = strcat('MA_cvAcc_','S',int2str(g),'.nii');
        H.descrip = sprintf('MA_cvLME_multi: model accuracy for session %s based on priors estimated from other sessions',int2str(g));
        spm_write_vol(H,reshape(S(g).Acc,m_dim));
        SPM.MACS.oosAcc(g) = H;
        H.fname   = strcat('MA_cvCom_','S',int2str(g),'.nii');
        H.descrip = sprintf('MA_cvLME_multi: model complexity for session %s based on priors estimated from other sessions',int2str(g));
        spm_write_vol(H,reshape(S(g).Com,m_dim));
        SPM.MACS.oosCom(g) = H;
    end;
end;

% Save SPM structure
%-------------------------------------------------------------------------%
save(strcat(SPM.swd,'/','SPM.mat'),'SPM');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);