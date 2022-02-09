function MS_SMM_BMS(BMS, mask, extent)
% _
% Determine Selected-Model Maps after Bayesian Model Selection
% FORMAT MS_SMM_BMS(BMS, mask, extent)
%     BMS    - a structure specifying an estimated BMS
%     mask   - a filepath to a binary masking image
%     extent - an integer, the voxel extent threshold
% 
% FORMAT MS_SMM_BMS(BMS, mask, extent) creates voxel-wise selected-
% model maps [2] for the models compared in BMS using a mask image
% indicating in-brain voxels and using an extent threshold defining the
% minimum number of neighboring voxels in order for a model to be selected
% in that cluster.
% 
% The output images written to disk are:
% - a continuous image for each model indicating
%   i)  whether it is selected (using non-NaN values);
%   ii) its likeliest frequency (where it is selected);
% - a selected model map indexing the model that was selected;
% - a map indicating the order position of the selected model.
% In case of family inference, these are selected-family maps.
% 
% Creation of selected-model maps is based on the likeliest frequencies
% from random-effects Bayesian model selection [1] and uses a cluster size
% correction algorithm that has an extent threshold (E) as input argument
% and works as follows:
% (1) Assign selected models based on posterior LFs.
% (2) For each voxel, determine the cluster size of the selected model.
% (3) In voxels where cluster size < E, select the next-best model.
% (4) Repeat (2) and (3) until all clusters have size >= E.
% 
% In this way, model clusters that are too small are iteratively removed.
% The default value for the input variable "mask" is the BMS mask image.
% The default value for the input variable "extent" is 10 voxels.
% 
% References:
% [1] Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ (2009):
%     "Bayesian model selection for group studies".
%     NeuroImage, vol. 46, pp. 1004-1017.
% [2] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469-489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/03/2015, 05:30 (V0.3/V10)
%  Last edit: 09/02/2022, 10:57 (V1.4/V20)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get BMS.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    BMS_mat = spm_select(1,'^BMS\.mat$','Select BMS.mat!');
    load(BMS_mat);
    MS_SMM_BMS(BMS);
    return
else
    BMS_dir = fileparts(BMS.fname);
    if nargin < 3 || isempty(extent), extent = 10; end;
    SMM_dir = strcat(BMS_dir,'/','MS_SMM_BMS_',num2str(extent),'/');
end;

% Set masking image if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(mask), mask = strcat(BMS_dir,'/',BMS.map.rfx.mask); end;

% Set extent threshold if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(extent), extent = 10; end;

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

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(BMS.map.rfx.alpha{1});          % alpha image header
V = prod(H.dim);                            % number of voxels

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_SMM_BMS: load');
spm_progress_bar('Init', 100, 'Load likeliest frequency maps...' , '');

% Load frequency maps
%-------------------------------------------------------------------------%
LF = zeros(M,V);
for i = 1:M
    lfm_hdr = spm_vol(BMS.map.rfx.lfm{i});
    lfm_img = spm_read_vols(lfm_hdr);
    lfm_img = reshape(lfm_img,[1 V]);
    LF(i,:) = lfm_img;
    spm_progress_bar('Set',(i/M)*100);
end;
clear lfm_hdr lfm_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Load mask image
%-------------------------------------------------------------------------%
m_hdr = spm_vol(mask);
m_img = spm_read_vols(m_hdr);
m_img = reshape(m_img,[1 V]);
m_ind = find(m_img~=0);


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Generate winning model maps
%-------------------------------------------------------------------------%
LF_max = repmat(max(LF),[M 1]); % maximal LF
LF_win = double(LF == LF_max);  % winning model

% Generate selected model map
%-------------------------------------------------------------------------%
M_sel = NaN(1,V);               % selected model
M_ord = NaN(1,V);               % selected model's
for j = 1:v                     % order position
    M_sel_ind       = find(LF_win(:,m_ind(j)));
    M_sel(m_ind(j)) = M_sel_ind(1);
    M_ord(m_ind(j)) = 1;
end;

% Generate sorted LF map
%-------------------------------------------------------------------------%
LF_sort = NaN(M,V);             % sorted LFs
for j = 1:v
    LF_sort(:,m_ind(j)) = sort(LF(:,m_ind(j)),'descend');
end;

% Generate cluster size map
%-------------------------------------------------------------------------%
M_sel_3D = reshape(M_sel,H.dim);
M_cls_3D = NaN(size(M_sel_3D));
for i = 1:M
    [cl_ind, num_cl] = spm_bwlabel(double(M_sel_3D==i),18);
    if num_cl > 0
        for l = 1:num_cl
            M_cls_3D(cl_ind==l) = sum(sum(sum(cl_ind==l)));
        end;
    end;
end;
M_cls   = reshape(M_cls_3D,[1 V]);
cls_min = min(M_cls(~isnan(M_cls)));

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_SMM_BMS: estimate');
spm_progress_bar('Init', 100, 'Eliminating sub-threshold clusters...', 'supra-threshold voxels');
spm_progress_bar('Set',(sum(M_cls>extent)/v)*100);

% Iteratively remove clusters
%-------------------------------------------------------------------------%
while cls_min < extent
    
    % Find small cluster voxels
    %---------------------------------------------------------------------%
    r_ind = find(M_cls<extent);
    num_r = length(r_ind);
    
    % Replace selected model there
    %---------------------------------------------------------------------%
    if num_r > 0
        for k = 1:num_r
            if M_ord(r_ind(k)) < M
                M_ord(r_ind(k)) = M_ord(r_ind(k)) + 1;
                M_sel_ind       = find(LF(:,r_ind(k))==LF_sort(M_ord(r_ind(k)),r_ind(k)));
                M_sel(r_ind(k)) = M_sel_ind(1);
            else
                M_ord(r_ind(k)) = M+1;      % Eliminate voxels that can not
                M_sel(r_ind(k)) = NaN;      % be settled, no matter what 
                m_ind(m_ind==r_ind(k)) = [];% model is placed there.
                v = v - 1;
            end;
        end;
    end;
    
    % Convert selected model to 3D
    %---------------------------------------------------------------------%
    M_sel_3D = reshape(M_sel,H.dim);
    M_cls_3D = NaN(size(M_sel_3D));
    
    % Label clusters in 3D image
    %---------------------------------------------------------------------%
    for i = 1:M
        [cl_ind, num_cl] = spm_bwlabel(double(M_sel_3D==i),18);
        if num_cl > 0
            for l = 1:num_cl
                M_cls_3D(cl_ind==l) = sum(sum(sum(cl_ind==l)));
            end;
        end;
    end;
    
    % Determine minimal cluster size
    %---------------------------------------------------------------------%
    M_cls   = reshape(M_cls_3D,[1 V]);
    cls_min = min(M_cls(~isnan(M_cls)));
    spm_progress_bar('Set',(sum(M_cls>extent)/v)*100);
    
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Generate binary model maps
%-------------------------------------------------------------------------%
M_bin = NaN(M,V);
for j = 1:v
    M_bin(M_sel(m_ind(j)),m_ind(j)) = 1;
end;

% Generate continuous model maps
%-------------------------------------------------------------------------%
M_con = M_bin .* LF;

% Sort models by voxels
%-------------------------------------------------------------------------%
mod_list = [[1:M]' zeros(M,1)];
for i = 1:M
    mod_list(i,2) = sum(M_sel==i);
end;
mod_list = sortrows(mod_list,-2);
sel_mods = sum(mod_list(:,2)>0);

% Retrieve model names
%-------------------------------------------------------------------------%
mod_names = cell(M,1);
for i = 1:M
    alpha_img = BMS.map.rfx.alpha{i};
    if ~isempty(strfind(alpha_img,'_model_'))
        mod_names{i} = alpha_img(1:strfind(alpha_img,'_model_')-1);
        MS_type      = 'model';
    elseif ~isempty(strfind(alpha_img,'_family_'))
        mod_names{i} = alpha_img(1:strfind(alpha_img,'_family_')-1);
        MS_type      = 'family';
    else
        mod_names{i} = strcat('GLM_',MF_int2str0(i,ceil(log10(M+1))));
        MS_type      = 'model';
    end;
end;
clear alpha_img


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_SMM_BMS: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(BMS.map.rfx.alpha{1});
if ~exist(SMM_dir,'dir')
    mkdir(SMM_dir);
end;

% Write images to disk
%-------------------------------------------------------------------------%
cd(SMM_dir);                                % selected model
H.fname   = strcat('MS_SMM_',MS_type(1:3),'.nii');
H.descrip = sprintf('MS_SMM_BMS: index of the selected %s', MS_type);
spm_write_vol(H,reshape(M_sel,H.dim));
H.fname   = strcat('MS_SMM_ord.nii');       % order position
H.descrip = sprintf('MS_SMM_BMS: position of the selected %s', MS_type);
spm_write_vol(H,reshape(M_ord,H.dim));
for i = 1:sel_mods                          % selected models
    H.fname   = strcat('MS_SMM_map','_pos_',MF_int2str0(i,ceil(log10(M+1))),'_',MS_type(1:3),'_',MF_int2str0(mod_list(i,1),ceil(log10(M+1))),'_',mod_names{mod_list(i,1)},'.nii');
    H.descrip = sprintf('MS_SMM_BMS: LF in voxels where %s %s has been selected', MS_type, mod_names{i});
    spm_write_vol(H,reshape(M_con(mod_list(i,1),:),H.dim));
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);