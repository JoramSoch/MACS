function B = MA_load_betas(SPM, m_ind)
% _
% Load Betas from General Linear Model
% FORMAT B = MA_load_betas(SPM, m_ind)
% 
%     SPM   - a structure specifying an estimated GLM
%     m_ind - a 1 x v vector indexing in-mask voxels
% 
%     B     - a p x v parameter matrix (p: regressors; v: voxels)
% 
% FORMAT B = MA_load_betas(SPM, m_ind) loads only in-mask time series
% belonging to the GLM specified by SPM and returns a parameter matrix.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 18:30 (V0.2/V6)
%  Last edit: 08/08/2017, 15:10 (V1.1/V17)


% Load mask if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(m_ind)
    [M m_dim m_ind] = MA_load_mask(SPM);
end;
clear M m_dim

% Get model dimensions
%-------------------------------------------------------------------------%
p = length(SPM.Vbeta);
v = numel(m_ind);
d = ceil(p/100);

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MA_load_betas: load');
spm_progress_bar('Init',100,'Load in-mask parameter estimates...','');

% Load parameter estimates
%-------------------------------------------------------------------------%
B = zeros(p,v);
for j = 1:p
    b_img  = spm_read_vols(SPM.Vbeta(j));
    B(j,:) = b_img(m_ind);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/p)*100); end;
end;
clear b_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');