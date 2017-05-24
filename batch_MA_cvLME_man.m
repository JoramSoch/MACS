function module = batch_MA_cvLME_man
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 16/03/2017, 21:55 (V0.99/V15)
%  Last edit: 17/03/2017, 06:25 (V0.99/V15)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MACS')); end


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Select SPM.mat
%-------------------------------------------------------------------------%
SPM_mat         = cfg_files;
SPM_mat.tag     = 'SPM_mat';
SPM_mat.name    = 'Select SPM.mat';
SPM_mat.help    = {'Select the SPM.mat file of a specified and/or estimated GLM.'};
SPM_mat.filter  = 'mat';
SPM_mat.ufilter = '^SPM\.mat$';
SPM_mat.num     = [1 1];

% Accuracy & Complexity
%-------------------------------------------------------------------------%
AnC        = cfg_menu;
AnC.tag    = 'AnC';
AnC.name   = 'Accuracy & Complexity';
AnC.help   = {'Choose whether model accuracy and model complexity are calculated in addition to the log model evidence.'};
AnC.labels = {'No', 'Yes'};
AnC.values = {0, 1};
AnC.val    = {0};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: cvLME (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_cvLME_man';
module.name = 'MA: calculate cvLME (manually)';
module.val  = {SPM_mat AnC};
module.help = {'Cross-Validated Log Model Evidence for General Linear Model'
               'Type "help MA_cvLME_multi" or "help MA_cvLME_single" for help.'};
module.prog = @run_module;
module.vout = @vout_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.SPM_mat{1});
AnC = logical(job.AnC);

% execute operation
MA_cvLME_multi(SPM,[],[],AnC);

% set output files
load(job.SPM_mat{1});
out.cvLME_nii = cellstr(strcat(SPM.swd,'/',SPM.MACS.cvLME.fname));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'cvLME map (voxel-wise image)';
dep(1).src_output = substruct('.','cvLME_nii');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});