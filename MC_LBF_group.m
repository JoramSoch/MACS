function MC_LBF_group(job)
% _
% Log Bayes Factors for General Linear Models
% FORMAT MC_LBF_group(job)
%     job    - SPM batch editor "MC: calculate LBF (manually)"
% 
% FORMAT MC_LBF_group(job) performs Bayesian model comparison (BMC) in a
% two-model model space and calculates the log Bayes factor (LBF), i.e.
% difference in log model evidence (LME), as well as the posterior
% probabilities (PP) [1,2] of the two models, for each subject
% as well as for the whole group of subjects [3].
% 
% Further information:
%     help ME_BMS_FFX
% 
% References:
% [1] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [2] Penny WD, Ridgway GR (2013): "Efficient Posterior
%     Probability Mapping Using Savage-Dickey Ratios". 
%     PLoS ONE, vol. 8, iss. 3, e59655.
% [3] Soch J, Allefeld C (2018): "MACS - a new SPM toolbox for model
%     assessment, comparison and selection". Journal of Neuroscience
%     Methods, in review. URL: https://www.biorxiv.org/content/early/2017/11/09/194365.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 17:05 (V0.99/V15)
%  Last edit: 09/03/2018, 11:55 (V1.2/V18)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get matlabbatch if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    design_mat = spm_select(1,'^*.mat','Select Batch Editor Job!');
    load(design_mat);
    job = matlabbatch{1}.spm.tools.MACS.MC_LBF_group_man;
    MS_LBF_group(job);
    return
else
    BMS_dir = job.dir{1};
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

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(job.models{1}{1}{1});%LME image header
V = prod(H.dim);                % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MC_LBF_group: load');
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');

% Load log model evidences
%-------------------------------------------------------------------------%
LME = zeros(M,V,N);             % M x V x N array of LMEs
for i = 1:N                     % subjects
    for j = 1:M                 % models
        lme_hdr = spm_vol(job.models{i}{j}{1});
        lme_img = spm_read_vols(lme_hdr);
        lme_img = reshape(lme_img,[1 V]);
        LME(j,:,i) = lme_img;
        spm_progress_bar('Set',(((i-1)*M+j)/(N*M))*100);
    end;
end;
clear lme_hdr lme_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
LME_1 = squeeze(LME(1,:,:));    % N x V matrix of LMEs for 1st model
if size(LME_1,2) == N, LME_1 = LME_1'; end;
[m_img m_hdr m_ind] = MS_create_mask(LME_1, H);
clear LME_1


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MC_LBF_group: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Calculate subject-level LBFs
%-------------------------------------------------------------------------%
LBF = zeros(2,V,N);
 PP = zeros(2,V,N);
for i = 1:N
    LBF(1,:,i) = LME(1,:,i) - LME(2,:,i);
    LBF(2,:,i) = LME(2,:,i) - LME(1,:,i);
    PP(1,:,i)  = exp(LBF(1,:,i)) ./ (exp(LBF(1,:,i)) + 1);
    PP(2,:,i)  = exp(LBF(2,:,i)) ./ (exp(LBF(2,:,i)) + 1);
end;

% Calculate group-level LBFs
%-------------------------------------------------------------------------%
LGBF = zeros(2,V);
 GPP = zeros(2,V);
LGBF(1,m_ind) = sum(LME(1,m_ind,:),3) - sum(LME(2,m_ind,:),3);
LGBF(2,m_ind) = sum(LME(2,m_ind,:),3) - sum(LME(1,m_ind,:),3);
 GPP(1,m_ind) = exp(LGBF(1,m_ind)) ./ (exp(LGBF(1,m_ind)) + 1);
 GPP(2,m_ind) = exp(LGBF(2,m_ind)) ./ (exp(LGBF(2,m_ind)) + 1);

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MC_LBF_group: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(job.models{1}{1}{1});

% Retrieve model names
%-------------------------------------------------------------------------%
mods = job.names';

% Write group-level images to disk
%-------------------------------------------------------------------------%
H.fname   = strcat('LBF_',mods{1},'_vs_',mods{2},'_group.nii');
H.descrip = sprintf('MC_LBF_group: log group Bayes factor (LGBF) in favor of %s',mods{1});
spm_write_vol(H,reshape(LGBF(1,:),H.dim));
H.fname   = strcat('LBF_',mods{2},'_vs_',mods{1},'_group.nii');
H.descrip = sprintf('MC_LBF_group: log group Bayes factor (LGBF) in favor of %s',mods{2});
spm_write_vol(H,reshape(LGBF(2,:),H.dim));
H.fname   = strcat('PP_',mods{1},'_group.nii');
H.descrip = sprintf('MC_LBF_group: posterior probability (PP) of %s',mods{1});
spm_write_vol(H,reshape(GPP(1,:),H.dim));
H.fname   = strcat('PP_',mods{2},'_group.nii');
H.descrip = sprintf('MC_LBF_group: posterior probability (PP) of %s',mods{2});
spm_write_vol(H,reshape(GPP(2,:),H.dim));

% Write subject-level images to disk
%-------------------------------------------------------------------------%
for i = 1:N
    H.fname   = strcat('LBF_',mods{1},'_vs_',mods{2},'_sub',MF_int2str0(i,ceil(log10(N+1))),'.nii');
    H.descrip = sprintf('MC_LBF_group: log Bayes factor (LBF) in favor of %s in subject %d',mods{1},i);
    spm_write_vol(H,reshape(LBF(1,:,i),H.dim));
    H.fname   = strcat('LBF_',mods{2},'_vs_',mods{1},'_sub',MF_int2str0(i,ceil(log10(N+1))),'.nii');
    H.descrip = sprintf('MC_LBF_group: log Bayes factor (LBF) in favor of %s in subject %d',mods{2},i);
    spm_write_vol(H,reshape(LBF(2,:,i),H.dim));
    H.fname   = strcat('PP_',mods{1},'_sub',MF_int2str0(i,ceil(log10(N+1))),'.nii');
    H.descrip = sprintf('MC_LBF_group: posterior probability (PP) of %s in subject %d',mods{1},i);
    spm_write_vol(H,reshape(PP(1,:,i),H.dim));
    H.fname   = strcat('PP_',mods{2},'_sub',MF_int2str0(i,ceil(log10(N+1))),'.nii');
    H.descrip = sprintf('MC_LBF_group: posterior probability (PP) of %s in subject %d',mods{2},i);
    spm_write_vol(H,reshape(PP(2,:,i),H.dim));
end;

% Write mask image to disk
%-------------------------------------------------------------------------%
spm_write_vol(m_hdr,reshape(m_img,H.dim));

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);