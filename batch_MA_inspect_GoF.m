function module = batch_MA_inspect_GoF
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 02:35 (V0.99/V15)
%  Last edit: 11/05/2017, 16:45 (V1.0/V16)


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
SPM_mat.help    = {'Select the SPM.mat file of an estimated GLM.'};
SPM_mat.filter  = 'mat';
SPM_mat.ufilter = '^SPM\.mat$';
SPM_mat.num     = [1 1];

% Type of plot
%-------------------------------------------------------------------------%
plot_type        = cfg_menu;
plot_type.tag    = 'plot_type';
plot_type.name   = 'Type of plot';
plot_type.help   = {'Choose between line graph (signals over time) and scatter plot (measured vs. predicted).'};
plot_type.labels = {'line graph', 'scatter plot'};
plot_type.values = {1, 2};
plot_type.val    = {1};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: inspect GoF
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_inspect_GoF';
module.name = 'MA: inspect goodness of fit';
module.val  = {SPM_mat plot_type};
module.help = {'Inspect Goodness of Fit in a General Linear Model'
               'Type "help MA_inspect_GoF" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.SPM_mat{1});
PlotTypes = {'linegraph','scatterplot'};
SPM.PlotType = PlotTypes{job.plot_type};

% execute operation
MA_inspect_GoF('Setup',SPM);

% set output files
out = [];