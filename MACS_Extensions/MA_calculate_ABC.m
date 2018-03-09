function MA_calculate_ABC(SPM)
% _
% Calculate ABC, A Bayesian Criterion

% Load design
%-------------------------------------------------------------------------%
X = SPM.xX.X;                   % design matrix
K = SPM.xX.K;                   % filtering matrix
W = SPM.xX.W;                   % whitening matrix
V = SPM.xVi.V;                  % non-sphericity

% Load data
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);
 Y              = MA_load_data(SPM,m_ind);

% Compute ABC
%-------------------------------------------------------------------------%
ABC        = NaN(size(M));
ABC(m_ind) = function_that_returns_voxel_wise_ABC(Y, X, K, W);

% Save ABC image
%-------------------------------------------------------------------------%
H         =  MA_init_header(SPM, false);
H.fname   = 'MA_ABC.nii';
H.descrip = 'MA_calculate ABC: A Bayesian Criterion';
spm_write_vol(H,reshape(ABC,m_dim));

% Save ABC header
%-------------------------------------------------------------------------%
SPM.MACS.ABC = H;
save(strcat(SPM.swd,'/','SPM.mat'),'SPM');