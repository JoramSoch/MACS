function [m_img m_hdr m_ind] = MS_create_mask(y_imgs, y_hdr)
% _
% Create Mask for Second-Level Analysis
% FORMAT [m_img m_hdr m_ind] = MS_create_mask(y_imgs, y_hdr)
% 
%     y_imgs - an N x V matrix of second-level data (N: number of subjects)
%     y_hdr  - a structure with exemplary header information for these data
% 
%     m_img  - a 1 x V mask vector (V: number of voxels)
%     m_hdr  - a structure with mask header information
%     m_ind  - a 1 x v vector indicating in-mask voxels
% 
% FORMAT [m_img m_hdr m_ind] = MS_create_mask(y_imgs, y_hdr) creates a mask
% image for the data y_imgs, such that all columns of this data matrix that
% have no NaN value are inside the mask, and a mask header that fits the
% input y_hdr. Then, it returns mask image and mask header along with
% indices of all in-mask voxels.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 09/12/2014, 13:40 (V0.2/V8)
%  Last edit: 17/08/2015, 23:20 (V0.3/V11)


% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_create_mask: save');

% Create mask image
%-------------------------------------------------------------------------%
m_img = double(~sign(sum(isnan(y_imgs),1)));
m_ind = find(m_img~=0);
m_hdr = y_hdr;

% Create mask header
%-------------------------------------------------------------------------%
m_hdr.fname   = 'mask.nii';
m_hdr.dt      = [spm_type('uint8') spm_platform('bigend')];
m_hdr.descrip = 'MS_create_mask: model selection mask';