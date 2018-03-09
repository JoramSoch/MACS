function module = batch_MS_DEF
% _
% Configure MATLAB Batch for MACS Toolbox

% Select MS.mat
%-------------------------------------------------------------------------%
MS_mat         = cfg_files;
MS_mat.tag     = 'MS_mat';
MS_mat.name    = 'Select MS.mat';
MS_mat.help    = {'Select the MS.mat file describing a model space of GLMs.'};
MS_mat.filter  = 'mat';
MS_mat.ufilter = '^MS\.mat$';
MS_mat.num     = [1 1];

% MS: DEF (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_DEF';
module.name = 'MS: perform DEF (automatic)';
module.val  = {MS_mat};
module.help = {'Discrete Evidence Fusion for General Linear Models'
               'Type "help MS_perform_DEF" for help.'};
module.prog = @run_module;
module.vout = @vout_module;

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% execute operation
load(job.MS_mat{1});
MS_perform_DEF(MS);
load(job.MS_mat{1});
out.DEF_mat = cellstr(strcat(MS.swd,'/','DEF.mat'));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'DEF results (DEF.mat file)';
dep(1).src_output = substruct('.','DEF_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});