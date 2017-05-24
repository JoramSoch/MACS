function MS_BMS_fams(BMS, mods, fams)
% _
% Bayesian Model Selection for General Linear Models (family inference)
% FORMAT MS_BMS_fams(BMS, mods, fams)
%     BMS  - a structure specifying an estimated BMS
%     mods - a 1 x M vector defining family affiliation (M: models)
%     fams - a 1 x F cell array defining family names (F: families)
% 
% FORMAT MS_BMS_fams(BMS, mods, fams) generates alpha maps for family
% inference. This procedure exploits the agglomerative property of the
% Dirichlet distribution [1] by summing Dirichlet-distributed random
% variables into new ones that are again Dirichlet-distributed with
% explicit parameter values [2].
% 
% The variable "mods" is a vector specifying for each model which family it
% belongs to. The variable "fams" is a cell array specifying the family
% names. It is required that length(fams) = max(mods).
% 
% Calling this function requires that the corresponding BMS has been
% estimated with a uniform prior over model families instead of a uniform
% prior over models. Otherwise, family inference will be biased [3].
% 
% NOTE: The MACS Toolbox discourages to use family-uniform priors (FUP) at
% group-level BMS and instead advocates to calculate log family evidences
% (LFE) from first-level LMEs [4]. This is implemented via the functions
% "MA_LFE_uniform" and "ME_MF_LFE".
% 
% Further information:
%     help MS_BMS_group
%     help MA_LFE_uniform
% 
% Exemplary usage:
%     MS_BMS_fams(BMS, [1 1 1 2 2], {'Fam_A' 'Fam_B'});
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, p. 1008, eq. 18.
% [2] Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB (2013):
%     "Bayesian Data Analysis". Chapman & Hall, 3rd edition, p. 583.
% [3] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709, p. 6.
% [4] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/11/2014, 14:45 (V0.2/V8)
%  Last edit: 17/03/2017, 09:45 (V0.99/V15)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get BMS directory
%-------------------------------------------------------------------------%
BMS_dir = fileparts(BMS.fname);

% Change to directory
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(BMS_dir);
catch
    BMS_dir = pwd;
end

% Get model parameters
%-------------------------------------------------------------------------%
N = length(BMS.map.rfx.data);               % number of subjects
M = length(BMS.map.rfx.alpha);              % number of models
F = length(fams); % = max(mods)             % number of families

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(BMS.map.rfx.alpha{1});          % alpha image header
V = prod(H.dim);                            % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_fams: load');
spm_progress_bar('Init', 100, 'Load model alpha maps...' , '');

% Load model alpha maps
%-------------------------------------------------------------------------%
Am = zeros(M,V);
for i = 1:M
    a_hdr   = spm_vol(BMS.map.rfx.alpha{i});
    a_img   = spm_read_vols(a_hdr);
    a_img   = reshape(a_img,[1 V]);
    Am(i,:) = a_img;                        % model-level parameters
    spm_progress_bar('Set',(i/M)*100);
end;
clear a_hdr a_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
[m_img m_hdr m_ind] = MS_create_mask(Am, H);


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_fams: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Compute family alpha maps
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Compute family alpha maps...', '');
Af = NaN(F,V);
for i = 1:F
    Af(i,:) = sum(Am(mods==i,:),1);         % family-level parameters
    spm_progress_bar('Set',(i/F)*100);
end;

% Compute likeliest frequencies
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Compute likeliest frequencies...', '');
LFM = NaN(F,V);
for j = 1:v
    LFM(:,m_ind(j)) = MD_Dir_mode(Af(:,m_ind(j)));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;

% Compute expected frequencies
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Compute expected frequencies...', '');
EFM = NaN(F,V);
for j = 1:v
    EFM(:,m_ind(j)) = MD_Dir_mean(Af(:,m_ind(j)));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;

% Compute exceedance probabilities
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Compute exceedance probabilities...', '');
EPM = NaN(F,V);
for j = 1:v
    EPM(:,m_ind(j)) = MD_Dir_exc_prob(Af(:,m_ind(j))')';
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMS_fams: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(BMS.map.rfx.alpha{1});

% Write images to disk
%-------------------------------------------------------------------------%
for i = 1:F
    H.fname   = strcat(fams{i},'_family_alpha.nii');
    H.descrip = 'MS_BMS_fams: family alpha maps';
    spm_write_vol(H,reshape(Af(i,:),H.dim));
    H.fname   = strcat(fams{i},'_family_LFM.nii');
    H.descrip = 'MS_BMS_fams: likeliest frequencies';
    spm_write_vol(H,reshape(LFM(i,:),H.dim));
    H.fname   = strcat(fams{i},'_family_EFM.nii');
    H.descrip = 'MS_BMS_fams: expected frequencies';
    spm_write_vol(H,reshape(EFM(i,:),H.dim));
    H.fname   = strcat(fams{i},'_family_EPM.nii');
    H.descrip = 'MS_BMS_fams: exceedance probabilities';
    spm_write_vol(H,reshape(EPM(i,:),H.dim));
end;

% Finalize BMS structure
%-------------------------------------------------------------------------%
BMS.map.rfx.family.mods = mods;
BMS.map.rfx.family.fams = fams;
save(strcat(BMS_dir,'/','BMS.mat'),'BMS');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);