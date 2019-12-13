function module = batch_MA_classic_ICs_man
% _
% Configure MATLAB Batch for "MA: classical ICs (manually)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 10:45 (V0.99/V15)
%  Last edit: 09/03/2018, 10:25 (V1.2/V18)


%=========================================================================%
% F I E L D S                                                             %
%=========================================================================%

% Select SPM.mat
%-------------------------------------------------------------------------%
SPM_mat         = cfg_files;
SPM_mat.tag     = 'SPM_mat';
SPM_mat.name    = 'Select SPM.mat';
SPM_mat.help    = {'Select the SPM.mat file of a specified and/or estimated GLM.'};
SPM_mat.filter  = 'mat';
SPM_mat.ufilter = '^SPM\.mat$';
SPM_mat.num     = [1 1];

% Information criterion
%--------------------------------------------------------------------------
ICs        = cfg_menu;
ICs.tag    = 'ICs';
ICs.name   = 'Information criterion';
ICs.help   = {'Select the classical information criterion that you want to calculate.'};
ICs.labels = {'Akaike information criterion (AIC)', ...
             'corrected Akaike information criterion (AICc)', ...
             'Bayesian information criterion (BIC)', ...
             'Deviance information criterion (DIC)', ...
             'Hannan-Quinn information criterion (HQC)', ...
             'Kullback information criterion (KIC)', ...
             'corrected Kullback information criterion (KICc)'};
ICs.values = {'AIC', 'AICc', 'BIC', 'DIC', 'HQC', 'KIC', 'KICc'};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: classic ICs (man)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_classic_ICs_man';
module.name = 'MA: classical ICs (manually)';
module.val  = {SPM_mat ICs};
module.help = {'Classical Information Criteria for General Linear Model'
               'Type "help MA_classic_ICs" for help.'};
module.prog = @run_module;
module.vout = @vout_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.SPM_mat{1});
ICs = {job.ICs};

% add working directory
if ~isfield(SPM,'swd')
    SPM.swd = fileparts(job.SPM_mat{1});
end;

% execute operation
MA_classic_ICs(SPM, [], ICs)

% set output files
load(job.SPM_mat{1});
eval(strcat('H = SPM.MACS.',ICs{1},';'));
out.IC_nii = cellstr(strcat(SPM.swd,'/',H.fname));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
ICs = {job.ICs};
dep(1)            = cfg_dep;
dep(1).sname      = sprintf('%s map (voxel-wise image)', ICs{1});
dep(1).src_output = substruct('.','IC_nii');
dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});