function MS_PPs_uniform(LMEs, mods, folder)
% _
% Posterior Model Probabilities assuming Uniform Prior over Models
% FORMAT MS_PPs_uniform(LMEs, mods, folder)
%     LMEs   - an M x 1 cell array specifying log model evidence maps
%     mods   - a  1 x M cell array indicating names of the models
%     folder - a string indicating the folder for the output images
% 
% FORMAT MS_PPs_uniform(LMEs, mods, folder) generates posterior probability
% maps [1] based on the log model evidence maps LMEs assuming uniform prior
% probabilities over models and saves them to sub-directories of the
% directory folder.
% 
% Further information:
%     help ME_MS_PPs
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/05/2018, 16:25 (V1.2/V18)
%  Last edit: 04/05/2018, 16:25 (V1.2/V18)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;
PPs_dir  = strcat(folder,'/','MS_PPs_uniform');

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});           % LME image header
M = numel(LMEs);                % number of models
V = prod(H.dim);                % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_PPs_uniform: load');
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');

% Load log model evidences
%-------------------------------------------------------------------------%
LME = zeros(M,V);               % M x V matrix of LMEs
for i = 1:M                     % subjects
    lme_hdr = spm_vol(LMEs{i});
    lme_img = spm_read_vols(lme_hdr);
    lme_img = reshape(lme_img,[1 V]);
    LME(i,:)= lme_img;
    spm_progress_bar('Set',(i/M)*100);
end;
clear lme_hdr lme_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_PPs_uniform: estimate');

% Calculate posterior probabilities
%-------------------------------------------------------------------------%
PPs = ME_MS_PPs(LME, (1/M)*ones(M,1));

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_PPs_uniform: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});
if ~exist(PPs_dir,'dir')
    mkdir(PPs_dir);
end;

% Write images to disk
%-------------------------------------------------------------------------%
cd(PPs_dir);
for i = 1:M
    H.fname   = strcat(mods{i},'_','PP.nii');
    H.descrip = 'MS_PPs_uniform: posterior probability assuming uniform prior';
    spm_write_vol(H,reshape(PPs(i,:),H.dim));
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);