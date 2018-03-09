function module = batch_MA_ABC
% _
% Configure MATLAB Batch for MACS Toolbox

% Select SPM.mat
%-------------------------------------------------------------------------%
SPM_mat         = cfg_files;
SPM_mat.tag     = 'SPM_mat';
SPM_mat.name    = 'Select SPM.mat';
SPM_mat.help    = {'Select the SPM.mat file of a specified and/or estimated GLM.'};
SPM_mat.filter  = 'mat';
SPM_mat.ufilter = '^SPM\.mat$';
SPM_mat.num     = [1 1];

% MA: ABC (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_ABC';
module.name = 'MA: calculate ABC (manually)';
module.val  = {SPM_mat};
module.help = {'A Bayesian Criterion for General Linear Model'
               'Type "help MA_calculate_ABC" for help.'};
module.prog = @run_module;
module.vout = @vout_module;

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% execute operation
load(job.SPM_mat{1});
MA_calculate_ABC(SPM);
load(job.SPM_mat{1});
out.ABC_nii = cellstr(strcat(SPM.swd,'/',SPM.MACS.ABC.fname));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'ABC map (voxel-wise image)';
dep(1).src_output = substruct('.','ABC_nii');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});