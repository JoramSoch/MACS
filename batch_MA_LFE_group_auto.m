function module = batch_MA_LFE_group_auto
% _
% Configure MATLAB Batch for "MA: calculate LFE (automatic)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 17:25 (V0.99/V15)
%  Last edit: 09/03/2018, 10:25 (V1.2/V18)


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


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: LFE group (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_LFE_group_auto';
module.name = 'MA: calculate LFE (automatic)';
module.val  = {MS_mat LME_map vector fam_names};
module.help = {'Log Family Evidence for a Family of General Linear Models'
               'Type "help MA_LFE_uniform" for help on this module.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.MS_mat{1});
mods = job.vector;
fams = job.names;

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
    MA_LFE_uniform(LMEs, mods, fams, folder);
end;

% set output files
out = [];