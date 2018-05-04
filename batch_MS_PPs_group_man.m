function module = batch_MS_PPs_group_man
% _
% Configure MATLAB Batch for "MS: calculate PPs (manually)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/05/2018, 16:35 (V1.2/V18)
%  Last edit: 04/05/2018, 16:35 (V1.2/V18)


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

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
subjects.num     = [1 Inf];

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
GLM_names.num    = [1 Inf];


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MS: PPs group (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_PPs_group_man';
module.name = 'MS: calculate PPs (manually)';
module.val  = {LME_maps GLM_names};
module.help = {'Posterior Probabilities within Model Space of General Linear Models'
               'Type "help MS_PPs_uniform" for help on this module.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
mods = job.names;

% execute operation
N = numel(job.models);
M = numel(job.models{1});
for i = 1:N
    LMEs = cell(M,1);
    for j = 1:M
        LMEs{j} = job.models{i}{j}{1};
    end;
    folder = strcat(fileparts(LMEs{1}),'/','..');
    MS_PPs_uniform(LMEs, mods, folder);
end;

% set output files
out = [];