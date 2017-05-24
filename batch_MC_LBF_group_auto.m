function module = batch_MC_LBF_group_auto
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 17:40 (V0.99/V15)
%  Last edit: 18/03/2017, 17:40 (V0.99/V15)


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


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MC: LBF group (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MC_LBF_group_man';
module.name = 'MC: calculate LBF (automatic)';
module.val  = {MS_mat LME_map};
module.help = {'Log Bayes Factors for General Linear Models'
               'Type "help MC_LBF_group" or "help ME_BMS_FFX" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.MS_mat{1});

% execute operation
[N,M] = size(MS.SPMs);
job.dir{1} = MS.swd;
job.names  = MS.GLMs;
job.models = cell(1,N);
for i = 1:N
    job.models{i} = cell(1,M);
    for j = 1:M
        load(MS.SPMs{i,j});
        eval(strcat('H = SPM.MACS.',job.LME_map,';'));
        job.models{i}{j} = cellstr(strcat(SPM.swd,'/',H.fname));
    end;
end;
MC_LBF_group(job);

% set output files
out = [];