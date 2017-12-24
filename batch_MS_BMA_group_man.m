function module = batch_MS_BMA_group_man
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 19:25 (V0.99/V15)
%  Last edit: 07/12/2017, 15:45 (V1.1/V17)


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

% Parameter matrix
%-------------------------------------------------------------------------%
para_mat         = cfg_entry;
para_mat.tag     = 'para_mat';
para_mat.name    = 'Parameter matrix';
para_mat.help    = {'The parameter matrix is an M x P matrix, where M is the number of models and P is the number of parameters.'
                    'Generally, the matrix indexes which regressor (column) in which GLM (row) belongs to which common mdoel parameter.'
                    'For example, in the case of three models (M = 3) and four parameters to be averaged (P = 4), it could look like this:'
                    '    1   2   3   4    '
                    '    1   3   5   7    '
                    '    1   4   7  10    '
                    'This would indicate that the 3rd/5th/7th regressor in the 1st/2nd/3rd model are the same model parameter.'
                    'The order of the models has to follow the order in which the cvLME maps are entered into the batch.'
                    'Note that, when a parameter is not included in a particular model and therefore estimated as 0 by this model,'
                    'this would be indicated by entering 0 into the parameter matrix. For example, it could look like this:'
                    '    1   2   3   4    '
                    '    1   3   5   7    '
                    '    1   4   7  10    '
                    '    1   4   0   0    '
                    'This would indicate that the 4th model does not have the 3rd/4th regressor from the 1st model.'
                    'Note that this variable can also be a 1 x P vector, if parameter indices are the same for all models.'};
para_mat.strtype = 'i';
para_mat.num     = [Inf Inf];

% Analysis name
%-------------------------------------------------------------------------%
analysis         = cfg_entry;
analysis.tag     = 'analysis';
analysis.name    = 'Analysis name';
analysis.help    = {'Enter a name without spaces to distinguish this analysis from others in the same model space.'
                    'If you are unsure what this means, just leave this variable at its default value.'};
analysis.strtype = 's';
analysis.num     = [1 Inf];
analysis.val     = {'BMA1'};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MS: BMA group (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_BMA_group_man';
module.name = 'MS: perform BMA (manually)';
module.val  = {LME_maps para_mat analysis};
module.help = {'Bayesian Model Averaging for General Linear Models'
               'Type "help MS_BMA_group" or "help MS_BMA_subject" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
params = job.para_mat;
method = strcat('ba','_',job.analysis);
job = rmfield(job,{'para_mat'});

% execute operation
MS_BMA_group(job,params,method);

% set output files
out = [];