function MA_LFE_uniform(LMEs, mods, fams, folder)
% _
% Log Family Evidence assuming Uniform Priors over and within Families
% FORMAT MA_LFE_uniform(LMEs, mods, fams, folder)
%     LMEs   - an M x 1 cell array specifying log model evidence maps
%     mods   - a  1 x M vector defining family affiliation (M: models)
%     fams   - a  1 x F cell array defining family names (F: families)
%     folder - a string indicating the folder for the output images
% 
% FORMAT MA_LFE_uniform(LMEs, mods, fams, folder) generates log family
% evidence maps based on the log model evidence maps LMEs assuming uniform
% prior probabilities over and within families and saves them to sub-
% directories of the directory folder.
% 
% The variable "mods" is a vector specifying for each model which family it
% belongs to. The variable "fams" is a cell array specifying the family
% names. It is required that length(fams) = max(mods).
% 
% Further information:
%     help ME_MF_LFE
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 05/02/2015, 12:20 (V0.3/V9)
%  Last edit: 07/12/2017, 20:50 (V1.1/V17)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;
LFE_dir  = strcat(folder,'/','MA_LFE_uniform');

% Get model parameters
%-------------------------------------------------------------------------%
M = length(LMEs);               % number of models
F = length(fams); % = max(mods) % number of families

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});           % LME image header
V = prod(H.dim);                % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_LFE_uniform: load');
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');

% Load log model evidences
%-------------------------------------------------------------------------%
LME = zeros(M,V);               % M x V array of LMEs
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
Finter = spm('FigName','MA_LFE_uniform: estimate');

% Calculate log family evidences
%-------------------------------------------------------------------------%
LFE = zeros(F,V);
for i = 1:F
    Mi = sum(mods==i);
    LFE(i,:) = ME_MF_LFE(LME(mods==i,:),(1/Mi)*ones(1,Mi));
end;
clear Mi

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_LFE_uniform: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});
if ~exist(LFE_dir,'dir')
    mkdir(LFE_dir);
end;

% Write images to disk
%-------------------------------------------------------------------------%
cd(LFE_dir);
for i = 1:F
    H.fname   = strcat(fams{i},'_','LFE.nii');
    H.descrip = 'MA_LFE_uniform: log family evidence assuming uniform priors';
    spm_write_vol(H,reshape(LFE(i,:),H.dim));
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);