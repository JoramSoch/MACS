function module = batch_MS_BMS_fams_auto
% _
% Configure MATLAB Batch for "MS: perform family BMS (automatic)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 09/03/2018, 10:25 (V1.2/V18)
%  Last edit: 09/02/2022, 11:24 (V1.4/V20)


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

% MS: BMS group (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MS_BMS_fams_auto';
module.name = 'MS: perform family BMS (automatic)';
module.val  = {MS_mat LME_map vector fam_names inf_meth EPs};
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
load(job.MS_mat{1});
method = job.inf_meth;
family = struct('mods',{job.vector},'fams',{job.names});
EPs    = logical(job.EPs);
job    = rmfield(job,{'vector','names','inf_meth','EPs'});

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
job = rmfield(job,{'MS_mat','LME_map'});
MS_BMS_group_fams(job,method,family,EPs);

% set output files
load(strcat(job.dir{1},'/','MS.mat'));
out.BMS_mat = cellstr(strcat(MS.swd,'/','BMS.mat'));

% Dependencies
%-------------------------------------------------------------------------%
function dep = vout_module(job)

% define dependencies
dep(1)            = cfg_dep;
dep(1).sname      = 'BMS results (BMS.mat file)';
dep(1).src_output = substruct('.','BMS_mat');
dep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});