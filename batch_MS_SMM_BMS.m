function module = batch_MS_SMM_BMS
% _
% Configure MATLAB Batch for MACS Toolbox
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 13:25 (V0.99/V15)
%  Last edit: 11/05/2017, 16:35 (V1.0/V16)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed, addpath(fullfile(spm('dir'),'toolbox','MACS')); end


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Select BMS.mat
%-------------------------------------------------------------------------%
BMS_mat         = cfg_files;
BMS_mat.tag     = 'BMS_mat';
BMS_mat.name    = 'Select BMS.mat';
BMS_mat.help    = {'Select the BMS.mat file of an estimated group-level BMS.'};
BMS_mat.filter  = 'mat';
BMS_mat.ufilter = '^BMS\.mat$';
BMS_mat.num     = [1 1];

% Extent threshold
%--------------------------------------------------------------------------
extent         = cfg_entry;
extent.tag     = 'extent';
extent.name    = 'Extent threshold';
extent.help    = {'Enter the voxel extent threshold, i.e. the minimum number of neighboring voxels in order for a model to be selected in that cluster.'};
extent.strtype = 'n';
extent.num     = [1 1];
extent.val     = {10};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MS: SMM BMS (opt)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_SMM_BMS';
module.name = 'MS: generate SMM from BMS';
module.val  = {BMS_mat extent};
module.help = {'Determine Selected-Model Maps after Bayesian Model Selection'
               'Type "help MS_SMM_BMS" for help on this module.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.BMS_mat{1});
extent = job.extent;

% execute operation
MS_SMM_BMS(BMS,[],extent);

% set output files
out = [];