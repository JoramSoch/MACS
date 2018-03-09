function module = batch_MF_visualize
% _
% Configure MATLAB Batch for "MF: visualize high-dimensional data"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 05:00 (V0.99/V15)
%  Last edit: 09/03/2018, 10:25 (V1.2/V18)


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Input images
%-------------------------------------------------------------------------%
data         = cfg_files;
data.tag     = 'data';
data.name    = 'Input images';
data.help    = {'Select your input images.'};
data.filter  = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];

% Overlay image
%-------------------------------------------------------------------------%
overlay         = cfg_files;
overlay.tag     = 'overlay';
overlay.name    = 'Overlay image';
overlay.help    = {'Select the overlay image.'};
overlay.filter  = 'image';
overlay.ufilter = '.*';
overlay.num     = [1 1];

% Overlay threshold
%-------------------------------------------------------------------------%
thresh         = cfg_entry;
thresh.tag     = 'thresh';
thresh.name    = 'Overlay threshold';
thresh.help    = {'Enter an overlay threshold, e.g. ''>0.5'' or ''==64''.'};
thresh.strtype = 's';
thresh.num     = [1 Inf];

% PlotType
%--------------------------------------------------------------------------
PlotType        = cfg_menu;
PlotType.tag    = 'PlotType';
PlotType.name   = 'PlotType';
PlotType.help   = {'Select the type of plot.'};
PlotType.labels = {'bar plot', 'line graph', 'matrix plot'};
PlotType.values = {'bar', 'plot', 'matrix'};
PlotType.val    = {'bar'};

% LineSpec
%-------------------------------------------------------------------------%
LineSpec         = cfg_entry;
LineSpec.tag     = 'LineSpec';
LineSpec.name    = 'LineSpec';
LineSpec.help    = {'Enter a LineSpec string.'};
LineSpec.strtype = 's';
LineSpec.num     = [1 Inf];
LineSpec.val     = {'b'};

% X-Axis Ticks
%-------------------------------------------------------------------------%
XTicks         = cfg_entry;
XTicks.tag     = 'XTicks';
XTicks.name    = 'X-Axis Ticks';
XTicks.help    = {'Enter X-Axis Ticks.'};
XTicks.strtype = 's';
XTicks.num     = [1 Inf];
XTicks.val     = {'{}'};

% Y-Axis Limits
%-------------------------------------------------------------------------%
YLimits         = cfg_entry;
YLimits.tag     = 'YLimits';
YLimits.name    = 'Y-Axis Limits';
YLimits.help    = {'Enter Y-Axis Limits.'};
YLimits.strtype = 's';
YLimits.num     = [1 Inf];
YLimits.val     = {'[]'};

% Title
%-------------------------------------------------------------------------%
Title         = cfg_entry;
Title.tag     = 'Title';
Title.name    = 'Title';
Title.help    = {'Enter the plot title.'};
Title.strtype = 's';
Title.num     = [1 Inf];
Title.val     = {'Title'};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MF: visualize
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MF_visualize';
module.name = 'MF: visualize high-dimensional data';
module.val  = {data overlay thresh PlotType LineSpec XTicks YLimits Title};
module.help = {'Visualization of High-Dimensional Data'
               'Type "help MF_visualize" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
xVis = job;
eval(strcat('xVis.XTicks = ',xVis.XTicks,';'));
xVis.YLimits = str2num(xVis.YLimits);

% execute operation
MF_visualize('Setup',xVis);

% set output files
out = [];