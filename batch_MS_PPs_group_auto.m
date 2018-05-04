function module = batch_MS_PPs_group_auto
% _
% Configure MATLAB Batch for "MS: calculate PPs (automatic)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/05/2018, 16:50 (V1.2/V18)
%  Last edit: 04/05/2018, 16:50 (V1.2/V18)


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

% MA: LFE group (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_PPs_group_auto';
module.name = 'MS: calculate PPs (automatic)';
module.val  = {MS_mat LME_map};
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
load(job.MS_mat{1});
mods = MS.GLMs;

% execute operation
[N,M] = size(MS.SPMs);
for i = 1:N
    LMEs = cell(M,1);
    for j = 1:M
        load(MS.SPMs{i,j});
        eval(strcat('H = SPM.MACS.',job.LME_map,';'));
        LMEs{j} = strcat(SPM.swd,'/',H.fname);
    end;
    folder = strcat(fileparts(LMEs{1}),'/','..');
    MS_PPs_uniform(LMEs, mods, folder);
end;

% set output files
out = [];