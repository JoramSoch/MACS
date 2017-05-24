function module = batch_MA_LFE_group_man
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 17:25 (V0.99/V15)
%  Last edit: 11/05/2017, 17:20 (V1.0/V16)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MACS')); end


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

% Family vector
%-------------------------------------------------------------------------%
vector         = cfg_entry;
vector.tag     = 'vector';
vector.name    = 'Family vector';
vector.help    = {'The family vector is a 1 x M vector defining family affiliation, where M is the number of models.'
                  'Generally, the number j in the i-th position indicates that the i-th model belongs to family j.'
                  'For example, the number 3 in the 5th position indicate that the 5th model belongs to family 3.'};
vector.strtype = 'n';
vector.num     = [1 Inf];

% Names
%-------------------------------------------------------------------------%
names         = cfg_entry;
names.tag     = 'names';
names.name    = 'Name';
names.help    = {'Enter the name of this model family of GLMs.'};
names.strtype = 's';
names.num     = [1 Inf];
names.val     = {'GLMs_'};

% Family names
%--------------------------------------------------------------------------
fam_names        = cfg_repeat;
fam_names.tag    = 'fam_names';
fam_names.name   = 'Family names';
fam_names.help   = {'The family names are a 1 x F cell array defining family names, where F is the number of families.'
                    'It is required that the numer of family names is equal to the highest number in the family vector,'
                    'i.e. it is required that length(names) = max(vector). Type "help MA_LFE_uniform" for further help.'};
fam_names.values = {names};
fam_names.num    = [1 Inf];


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: LFE group (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_LFE_group_man';
module.name = 'MA: calculate LFE (manually)';
module.val  = {LME_maps vector fam_names};
module.help = {'Log Family Evidence for a Family of General Linear Models'
               'Type "help MA_LFE_uniform" for help on this module.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
mods = job.vector;
fams = job.names;

% execute operation
N = numel(job.models);
M = numel(job.models{1});
for i = 1:N
    LMEs = cell(M,1);
    for j = 1:M
        LMEs{j} = job.models{i}{j}{1};
    end;
    folder = strcat(fileparts(LMEs{1}),'/','..');
    MA_LFE_uniform(LMEs, mods, fams, folder);
end;

% set output files
out = [];