function [M m_dim m_ind] = MA_load_mask(SPM)
% _
% Load Mask from General Linear Model
% FORMAT [M m_dim m_ind] = MA_load_mask(SPM)
% 
%     SPM   - a structure specifying an estimated GLM
% 
%     M     - a 1 x V mask vector (V: number of voxels)
%     m_dim - a 1 x 3 vector indicating mask dimensions
%     m_ind - a 1 x v vector indicating in-mask voxels
% 
% FORMAT [M m_dim m_ind] = MA_load_mask(SPM) loads the mask image belonging
% to the GLM specified by SPM and returns the mask, its dimensions and
% indices of all in-mask voxels.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 18:20 (V0.2/V6)
%  Last edit: 27/11/2014, 17:00 (V0.2/V8)


% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_load_mask: load');

% Load mask image
%-------------------------------------------------------------------------%
m_dim = SPM.VM.dim;
m_img = spm_read_vols(SPM.VM);
m_img = reshape(m_img,[1 prod(m_dim)]);
m_ind = find(m_img~=0);
M     = m_img;