function MS_BMS_group(job, method, family, EPs)
% _
% Bayesian Model Selection for General Linear Models (model inference)
% FORMAT MS_BMS_group(job, method, family, EPs)
%     job    - SPM batch editor "MS: perform BMS (manually)"
%     method - a string indicating which methods to use
%     family - a structure with the following fields:
%     o mods - a 1 x M vector defining family affiliation (M: models)
%     o fams - a 1 x F cell array defining family names (F: families)
%     EPs    - a logical indicating exceedance probability calculation
% 
% FORMAT MS_BMS_group(job, method, family, EPs) generates Bayesian model
% selection maps according to the batch editor job and using the method
% indicated by method with family inference indicated by family.
% 
% The input variable "method" is a string indicating which method to use:
%     If method is 'FFX',    then a fixed effects model is estimated.
%     If method is 'RFX-VB', then a Variational Bayes approach is taken.
%     If method is 'RFX-GS', then a Gibbs Sampling approach is taken.
% 
% The default for this variable is 'RFX-VB' which means that the second-
% level model over log model evidence is estimated using Variational Bayes
% [1,2]. This is recommended if a large space containing lots of voxels is
% analyzed. Setting the method to 'RFX-GS' invokes estimation via Gibbs
% Sampling [3] which is more time-consuming, but also more precise.
% 
% The input variable "family" is an optional structure that specifies
% family inference. The field "mods" is a vector specifying for each model
% which family it belongs to. The field "fams" is a cell array specifying
% the family names. It is required that length(fams) = max(mods). By
% default, this variable is empty.
% 
% When this variable is empty (as by default), a uniform prior over models
% is applied. If this variable is non-empty, this invokes a family-wise
% Bayesian model selection using the function "MS_BMS_group_fams". This
% means that first, log family evidences (LFE) are calculated from first-
% level LMEs [4], and then, Bayesian model selection with a uniform prior
% over families is applied.
% 
% The input variable "EPs" is a logical indicating whether exceedance
% probabilities (EP) are calculated and written to images in case method is
% 'RFX-VB'. If the number of models is very large, it is advisable to only
% calculate EPs at the family level. By default, this variable is false.
% 
% Further information:
%     help ME_BMS_FFX
%     help ME_BMS_RFX_VB
%     help ME_BMS_RFX_GS
%     help MD_Dir_exc_prob
% 
% Exemplary usage:
%     MS_BMS_group(matlabbatch, 'RFX-VB', [], true);
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Rosa MJ, Bestmann S, Harrison L, Penny W (2010):
%     "Bayesian model selection maps for group studies".
%     NeuroImage, vol. 49, pp. 217-224.
% [3] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [4] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469�489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 05/12/2014, 12:15 (V0.2/V8)
%  Last edit: 09/02/2022, 10:21 (V1.4/V20)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get matlabbatch if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    design_mat = spm_select(1,'^*.mat','Select Batch Editor Job!');
    load(design_mat);
    job = matlabbatch{1}.spm.tools.MACS.MS_BMS_group_man;
    MS_BMS_group(job);
    return
else
    BMS_dir = job.dir{1};
end;

% Set inference method if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(method), method = 'RFX-VB'; end;

% Set family inference if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(family), family = []; end;

% Inactivate EPs if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(EPs), EPs = false; end;

% Call other function if families
%-------------------------------------------------------------------------%
if ~isempty(family)
    MS_BMS_group_fams(job, method, family, EPs);
    return
end;

% Change to directory if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(BMS_dir);
catch
    BMS_dir = strcat(pwd,'/');
end

% Get model parameters
%-------------------------------------------------------------------------%
N = numel(job.models);          % number of subjects
M = numel(job.names);           % number of models
S = 1; % already eliminated     % number of sessions

% Get family dimensions
%-------------------------------------------------------------------------%
if ~isempty(family)
    F  = length(family.fams);
    Mk = sum(repmat([1:F]',[1 M])==repmat(family.mods,[F 1]),2)';
    Mj = Mk(family.mods);
  % F  - the number of families [1 x 1]
  % Mk - the number of models in family k [1 x F]
  % Mj - the number of models in family that model j belongs to [1 x M]
end;

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(job.models{1}{1}{1});           % LME image header
V = prod(H.dim);                            % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group: load');
spm_progress_bar('Init', 100, 'Load log model evidences...', '');

% Load log model evidences
%-------------------------------------------------------------------------%
LME = zeros(N,M,V);             % N x M x V array of LMEs
for i = 1:N                     % subjects
    for j = 1:M                 % models
        lme_hdr = spm_vol(job.models{i}{j}{S});
        lme_img = spm_read_vols(lme_hdr);
        lme_img = reshape(lme_img,[1 1 V]);
        LME(i,j,:) = lme_img;
        % NOTE: This group-level BMS operates on subject- and model-, but
        % not session-wise LMEs, because session-wise (oos)LMEs have
        % already been summed up during first-level model assessment.
        spm_progress_bar('Set',(((i-1)*M+j)/(N*M))*100);
    end;
end;
clear lme_hdr lme_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
LMEs = reshape(LME,[N*M, V]);   % (N*M) x V matrix of LMEs
[m_img m_hdr m_ind] = MS_create_mask(LMEs, H);
clear LMEs


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Fixed Effects Inference (Bayes' Rule)
%-------------------------------------------------------------------------%
if strcmp(method,'FFX')
    % prior and posterior
    if isempty(family)
        prior = 1/M * ones(1,M);
    else
        prior = 1/F * 1./Mj;
    end;
    PPM   = NaN(M,V);           % posterior probability maps
    % voxel-wise estimation
    spm_progress_bar('Init', 100, 'Estimate posterior probabilities...', '');
    for j = 1:v
        [LGBF, post] = ME_BMS_FFX(LME(:,:,m_ind(j)), prior);
        PPM(:,m_ind(j)) = post';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
    clear LGBF
end;

% Random Effects Inference (Variational Bayes)
%-------------------------------------------------------------------------%
if strcmp(method,'RFX-VB')
    % prior and posterior
    if isempty(family)
        prior = ones(1,M);
    else
        prior = 1./Mj;
    end;
    alpha = NaN(M,V);           % alpha parameter maps
    LFM   = NaN(M,V);           % likeliest frequency maps
    EFM   = NaN(M,V);           % expected frequency maps
    EPM   = NaN(M,V);           % exceedance probability maps
    % voxel-wise estimation
    spm_progress_bar('Init', 100, 'Estimate posterior distribution...', '');
    for j = 1:v
        alpha(:,m_ind(j)) = ME_BMS_RFX_VB(LME(:,:,m_ind(j)), prior)';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
    % voxel-wise computation
    spm_progress_bar('Init', 100, 'Compute posterior frequencies...', '');
    for j = 1:v
        LFM(:,m_ind(j)) = MD_Dir_mode(alpha(:,m_ind(j)));
        EFM(:,m_ind(j)) = MD_Dir_mean(alpha(:,m_ind(j)));
        if EPs, EPM(:,m_ind(j)) = MD_Dir_exc_prob(alpha(:,m_ind(j))')'; end;
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;

% Random Effects Inference (Gibbs Sampling)
%-------------------------------------------------------------------------%
if strcmp(method,'RFX-GS')
    % prior and posterior
    if isempty(family)
        prior = ones(1,M);
    else
        prior = 1./Mj;
    end;
    alpha = NaN(M,V);           % alpha parameter maps
    EFM   = NaN(M,V);           % expected frequency maps
    EPM   = NaN(M,V);           % exceedance probabilties maps
    % voxel-wise sampling
    spm_progress_bar('Init', 100, 'Sample posterior distribution...', '');
    for j = 1:v
        [alpha_post, exp_freq, exc_prob] = ME_BMS_RFX_GS(LME(:,:,m_ind(j)), prior);
        alpha(:,m_ind(j)) = alpha_post';
        EFM(:,m_ind(j))   = exp_prob';
        EPM(:,m_ind(j))   = exc_prob';
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_group: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(job.models{1}{1}{1});

% Retrieve model names
%-------------------------------------------------------------------------%
mods = job.names';

% Write images to disk
%-------------------------------------------------------------------------%
for i = 1:M
    if strcmp(method,'FFX')     % Bayes' Theorem
        H.fname   = strcat(mods{i},'_model_PPM.nii');
        H.descrip = 'MS_BMS_group: posterior probability maps';
        spm_write_vol(H,reshape(PPM(i,:),H.dim));
        BMS.map.ffx.ppm{i} = H.fname;
    end;
    if strcmp(method,'RFX-VB')  % Variational Bayes
        H.fname   = strcat(mods{i},'_model_alpha.nii');
        H.descrip = 'MS_BMS_group: alpha parameters';
        spm_write_vol(H,reshape(alpha(i,:),H.dim));
        BMS.map.rfx.alpha{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_LFM.nii');
        H.descrip = 'MS_BMS_group: likeliest frequencies';
        spm_write_vol(H,reshape(LFM(i,:),H.dim));
        BMS.map.rfx.lfm{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EFM.nii');
        H.descrip = 'MS_BMS_group: expected frequencies';
        spm_write_vol(H,reshape(EFM(i,:),H.dim));
        BMS.map.rfx.efm{i} = H.fname;
        if EPs
            H.fname   = strcat(mods{i},'_model_EPM.nii');
            H.descrip = 'MS_BMS_group: exceedance probabilities';
            spm_write_vol(H,reshape(EPM(i,:),H.dim));
            BMS.map.rfx.epm{i} = H.fname;
        end;
    end;
    if strcmp(method,'RFX-GS')  % Gibbs Sampling
        H.fname   = strcat(mods{i},'_model_alpha.nii');
        H.descrip = 'MS_BMS_group: alpha parameters';
        spm_write_vol(H,reshape(alpha(i,:),H.dim));
        BMS.map.rfx.alpha{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EFM.nii');
        H.descrip = 'MS_BMS_group: expected frequencies';
        spm_write_vol(H,reshape(EFM(i,:),H.dim));
        BMS.map.rfx.efm{i} = H.fname;
        H.fname   = strcat(mods{i},'_model_EPM.nii');
        H.descrip = 'MS_BMS_group: exceedance probabilities';
        spm_write_vol(H,reshape(EPM(i,:),H.dim));
        BMS.map.rfx.epm{i} = H.fname;
    end;
end;

% Write mask image to disk
%-------------------------------------------------------------------------%
spm_write_vol(m_hdr,reshape(m_img,H.dim));

% Finalize BMS structure
%-------------------------------------------------------------------------%
if strcmp(method,'FFX')
    BMS.map.ffx.data  = job.models;
    BMS.map.ffx.mask  = m_hdr.fname;
    BMS.map.ffx.prior = prior;
end;
if strcmp(method,'RFX-VB') || strcmp(method,'RFX-GS')
    BMS.map.rfx.data  = job.models;
    BMS.map.rfx.mask  = m_hdr.fname;
    BMS.map.rfx.prior = prior;
end;
BMS.fname = strcat(BMS_dir,'/','BMS.mat');
save(strcat(BMS_dir,'/','BMS.mat'),'BMS');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);