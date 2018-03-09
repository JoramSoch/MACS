function module = batch_MC_LBF_group_man
% _
% Configure MATLAB Batch for "MC: calculate LBF (manually)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 17:31 (V0.99/V15)
%  Last edit: 09/03/2018, 10:25 (V1.2/V18)


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Directory
%-------------------------------------------------------------------------%
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select the directory into which the LBF images will be saved.'};
dir.filter  = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% Models
%--------------------------------------------------------------------------
models         = cfg_files;
models.tag     = 'models';
models.name    = 'Model';
models.help    = {'Select the log model evidence map for this GLM in this subject.'
                  'Models should be specfied in the same order for each subject.'};
models.filter  = 'image';
models.ufilter = '.*';
models.num     = [1 1];

% Subjects
%--------------------------------------------------------------------------
subjects         = cfg_repeat;
subjects.tag     = 'subjects';
subjects.name    = 'Subject';
subjects.help    = {'Select the log model evidence maps of the GLMs for this subject.'
                    'Models should be specfied in the same order for each subject.'};
subjects.values  = {models};
subjects.num     = [2 2];

% Select LME maps
%-------------------------------------------------------------------------%
LME_maps         = cfg_repeat;
LME_maps.tag     = 'LME_maps';
LME_maps.name    = 'Select LME maps';
LME_maps.help    = {'Select the log model evidence map for each subject and model.'};
LME_maps.values  = {subjects};
LME_maps.num     = [1 Inf];

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
GLM_names.help   = {'Enter a name for each model.'}';
GLM_names.values = {names};
GLM_names.num    = [2 2];


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MC: LBF group (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MC_LBF_group_man';
module.name = 'MC: calculate LBF (manually)';
module.val  = {dir LME_maps GLM_names};
module.help = {'Log Bayes Factors for General Linear Models'
               'Type "help MC_LBF_group" or "help ME_BMS_FFX" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% execute operation
MC_LBF_group(job);

% set output files
out = [];