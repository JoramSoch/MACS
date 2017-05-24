function module = batch_MA_model_space
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 01:35 (V0.99/V15)
%  Last edit: 17/03/2017, 02:20 (V0.99/V15)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MACS')); end


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Directory
%-------------------------------------------------------------------------%
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select the directory into which the MS.mat file will be saved.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% Models
%--------------------------------------------------------------------------
models         = cfg_files;
models.tag     = 'models';
models.name    = 'Model';
models.help    = {'Select the SPM.mat file of this GLM in this subject.'
                  'Models should be specfied in the same order for each subject.'};
models.filter  = 'mat';
models.ufilter = '^SPM\.mat$';
models.num     = [1 1];

% Subjects
%--------------------------------------------------------------------------
subjects         = cfg_repeat;
subjects.tag     = 'subjects';
subjects.name    = 'Subject';
subjects.help    = {'Select the SPM.mat files of the GLMs for this subject.'
                    'Models should be specfied in the same order for each subject.'};
subjects.values  = {models};
subjects.num     = [1 Inf];

% Select SPM.mats
%-------------------------------------------------------------------------%
SPM_mats         = cfg_repeat;
SPM_mats.tag     = 'SPM_mats';
SPM_mats.name    = 'Select SPM.mats';
SPM_mats.help    = {'Select the SPM.mat file for each subject and model.'};
SPM_mats.values  = {subjects};
SPM_mats.num     = [1 Inf];

% Names
%-------------------------------------------------------------------------%
names         = cfg_entry;
names.tag     = 'names';
names.name    = 'Name';
names.help    = {'Enter the name of this GLM.'};
names.strtype = 's';
names.num     = [1 Inf];
names.val     = {'GLM_'};

% Enter GLM names
%--------------------------------------------------------------------------
GLM_names        = cfg_repeat;
GLM_names.tag    = 'GLM_names';
GLM_names.name   = 'Enter GLM names';
GLM_names.help   = {'Enter a name for each model.'};
GLM_names.values = {names};
GLM_names.num    = [1 Inf];


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: model space
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_model_space';
module.name = 'MA: define model space';
module.val  = {dir SPM_mats GLM_names};
module.help = {'Model Space of General Linear Models'
               'Create MS.mat file specifying a model space of GLMs.'};
module.prog = @run_module;
module.vout = @vout_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% execute operation
MA_model_space(job);

% set output files
out.MS_mat = cellstr(strcat(job.dir{1},'/','MS.mat'));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'model space (MS.mat file)';
dep(1).src_output = substruct('.','MS_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});