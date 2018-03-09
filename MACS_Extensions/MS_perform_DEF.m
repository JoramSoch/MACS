function MS_perform_DEF(MS)
% _
% Perform DEF, Discrete Evidence Fusion

% Get model parameters
%-------------------------------------------------------------------------%
N = size(MS.SPMs,1);            % number of subjects
M = size(MS.SPMs,2);            % number of models

% Get image dimensions
%-------------------------------------------------------------------------%
load(MS.SPMs{1,1});             % first subject/model
H = spm_vol(SPM.MACS.ABC);      % ABC image header
V = prod(H.dim);                % number of voxels

% Load ABC images
%-------------------------------------------------------------------------%
ABC = zeros(N,M,V);             % N x M x V array of ABCs
for i = 1:N                     % subjects
    for j = 1:M                 % models
        load(MS.SPMs{i,j});
        abc_hdr    = spm_vol(SPM.MACS.ABC);
        abc_img    = spm_read_vols(abc_hdr);
        ABC(i,j,:) = reshape(abc_img,[1 1 V]);
    end;
end;

% Create mask image
%-------------------------------------------------------------------------%
ABC_1 = squeeze(ABC(:,1,:));    % N x V matrix of ABCs for 1st model
[m_img m_hdr m_ind] = MS_create_mask(ABC_1, H);
 v    = numel(m_ind);

% Perform DEF analysis
%-------------------------------------------------------------------------%
DEF = NaN(M,V);
for j = 1:v
    DEF(:,m_ind(j)) = function_that_performs_voxel_wise_DEF(ABC(:,:,j))';
end;

% Save DEF maps
%-------------------------------------------------------------------------%
cd(MS.swd);
for i = 1:M
    H.fname   = strcat(MS.GLMs{i},'_DEF.nii');
    H.descrip = 'MS_perform DEF: Discrete Evidence Fusion';
    spm_write_vol(H,reshape(DEF(i,:),H.dim));
end;