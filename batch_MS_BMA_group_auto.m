function module = batch_MS_BMA_group_auto
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 22:20 (V0.99/V15)
%  Last edit: 11/05/2017, 19:05 (V1.0/V16)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MACS')); end


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Select MS.mat
%-------------------------------------------------------------------------%
MS_mat         = cfg_files;
MS_mat.tag     = 'MS_mat';
MS_mat.name    = 'Select MS.mat';
MS_mat.help    = {'Select the MS.mat file describing a model space of GLMs.'};
MS_mat.filter  = 'mat';
MS_mat.ufilter = '^MS\.mat$';
MS_mat.num     = [1 1];

% Enter LME map
%-------------------------------------------------------------------------%
LME_map         = cfg_entry;
LME_map.tag     = 'LME_map';
LME_map.name    = 'Enter LME map';
LME_map.help    = {'Enter the name of the log model evidence map that you want to use.'
                   'In the SPM.mat, LME maps will then be accessed via "SPM.MACS.[input].fname".'
                   'If you are unsure what this means, just leave this variable at its default value.'};
LME_map.strtype = 's';
LME_map.num     = [1 Inf];
LME_map.val     = {'cvLME'};

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
para_mat.strtype = 'n';
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

% MS: BMA group (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_BMA_group_auto';
module.name = 'MS: perform BMA (automatic)';
module.val  = {MS_mat LME_map para_mat analysis};
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
load(job.MS_mat{1});
params = job.para_mat;
method = strcat('ba','_',job.analysis);
job = rmfield(job,{'para_mat'});

% execute operation
[N,M] = size(MS.SPMs);
job.models = cell(1,N);
for i = 1:N
    job.models{i} = cell(1,M);
    for j = 1:M
        load(MS.SPMs{i,j});
        eval(strcat('H = SPM.MACS.',job.LME_map,';'));
        job.models{i}{j} = cellstr(strcat(SPM.swd,'/',H.fname));
    end;
end;
job = rmfield(job,{'MS_mat','LME_map'});
MS_BMA_group(job,params,method);

% set output files
out = [];