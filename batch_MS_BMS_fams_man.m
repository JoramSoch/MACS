function module = batch_MS_BMS_fams_man
% _
% Configure MATLAB Batch for "MS: perform family BMS (manually)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 09/03/2018, 10:25 (V1.2/V18)
%  Last edit: 09/02/2022, 11:37 (V1.4/V20)


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Directory
%-------------------------------------------------------------------------%
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select the directory into which the BMS.mat file will be saved.'};
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
                  'For example, the number 3 in the 5th position would indicate that the 5th model belongs to family 3.'};
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

% Inference method
%--------------------------------------------------------------------------
inf_meth        = cfg_menu;
inf_meth.tag    = 'inf_meth';
inf_meth.name   = 'Inference method';
inf_meth.help   = {'Select inference method using which Bayesian model selection should be performed.'};
inf_meth.labels = {'Fixed Effects (FFX)', ...
                   'Random Effects with Variational Bayes (RFX-VB)', ...
                   'Random Effects with Gibbs Sampling (RFX-GS)'};
inf_meth.values = {'FFX', 'RFX-VB', 'RFX-GS'};
inf_meth.val    = {'RFX-VB'};

% Exceedance probabilities
%-------------------------------------------------------------------------%
EPs        = cfg_menu;
EPs.tag    = 'EPs';
EPs.name   = 'Exceedance probabilities';
EPs.help   = {'Choose whether exceedance probabilities are calculated in addition to expected and likeliest frequencies.'};
EPs.labels = {'No', 'Yes'};
EPs.values = {0, 1};
EPs.val    = {0};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MS: BMS group (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_BMS_fams_man';
module.name = 'MS: perform family BMS (manually)';
module.val  = {dir LME_maps vector fam_names inf_meth EPs};
module.help = {'Family-Wise Bayesian Model Selection for General Linear Models'
               'Type "help MS_BMS_group_fams" or "help ME_BMS_RFX_VB" for help.'};
module.prog = @run_module;
module.vout = @vout_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
method = job.inf_meth;
family = struct('mods',{job.vector},'fams',{job.names});
EPs    = logical(job.EPs);
job    = rmfield(job,{'vector','names','inf_meth','EPs'});

% execute operation
job.names = cell(1,numel(job.models{1}));
MS_BMS_group_fams(job,method,family,EPs);

% set output files
out.BMS_mat = cellstr(strcat(job.dir{1},'/','BMS.mat'));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'BMS results (BMS.mat file)';
dep(1).src_output = substruct('.','BMS_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});