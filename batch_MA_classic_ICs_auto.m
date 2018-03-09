function module = batch_MA_classic_ICs_auto
% _
% Configure MATLAB Batch for "MA: classical ICs (automatic)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/03/2017, 11:45 (V0.99/V15)
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

% AIC
%-------------------------------------------------------------------------%
AIC        = cfg_menu;
AIC.tag    = 'AIC';
AIC.name   = 'AIC';
AIC.help   = {'Akaike information criterion (AIC)'};
AIC.labels = {'No', 'Yes'};
AIC.values = {0, 1};
AIC.val    = {1};

% AICc
%-------------------------------------------------------------------------%
AICc        = cfg_menu;
AICc.tag    = 'AICc';
AICc.name   = 'AICc';
AICc.help   = {'corrected Akaike information criterion (AICc)'};
AICc.labels = {'No', 'Yes'};
AICc.values = {0, 1};
AICc.val    = {0};

% BIC
%-------------------------------------------------------------------------%
BIC        = cfg_menu;
BIC.tag    = 'BIC';
BIC.name   = 'BIC';
BIC.help   = {'Bayesian information criterion (BIC)'};
BIC.labels = {'No', 'Yes'};
BIC.values = {0, 1};
BIC.val    = {1};

% DIC
%-------------------------------------------------------------------------%
DIC        = cfg_menu;
DIC.tag    = 'DIC';
DIC.name   = 'DIC';
DIC.help   = {'Deviance information criterion (DIC)'};
DIC.labels = {'No', 'Yes'};
DIC.values = {0, 1};
DIC.val    = {0};

% HQC
%-------------------------------------------------------------------------%
HQC        = cfg_menu;
HQC.tag    = 'HQC';
HQC.name   = 'HQC';
HQC.help   = {'Hannan-Quinn information criterion (HQC)'};
HQC.labels = {'No', 'Yes'};
HQC.values = {0, 1};
HQC.val    = {0};

% KIC
%-------------------------------------------------------------------------%
KIC        = cfg_menu;
KIC.tag    = 'KIC';
KIC.name   = 'KIC';
KIC.help   = {'Kullback information criterion (KIC)'};
KIC.labels = {'No', 'Yes'};
KIC.values = {0, 1};
KIC.val    = {0};

% KICc
%-------------------------------------------------------------------------%
KICc        = cfg_menu;
KICc.tag    = 'KICc';
KICc.name   = 'KICc';
KICc.help   = {'corrected Kullback information criterion (KICc)'};
KICc.labels = {'No', 'Yes'};
KICc.values = {0, 1};
KICc.val    = {0};

% Information criteria
%-------------------------------------------------------------------------%
ICs      = cfg_branch;
ICs.tag  = 'ICs';
ICs.name = 'Information criteria';
ICs.val  = {AIC AICc BIC DIC HQC KIC KICc};
ICs.help = {'Select the classical information criteria that you want to calculate.'};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: classic ICs (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_classic_ICs_auto';
module.name = 'MA: classical ICs (automatic)';
module.val  = {MS_mat ICs};
module.help = {'Classical Information Criteria for General Linear Model'
               'Type "help MA_classic_ICs" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.MS_mat{1});
ICs = [];
if job.ICs.AIC  == 1, ICs = [ICs {'AIC'}];  end;
if job.ICs.AICc == 1, ICs = [ICs {'AICc'}]; end;
if job.ICs.BIC  == 1, ICs = [ICs {'BIC'}];  end;
if job.ICs.DIC  == 1, ICs = [ICs {'DIC'}];  end;
if job.ICs.HQC  == 1, ICs = [ICs {'HQC'}];  end;
if job.ICs.KIC  == 1, ICs = [ICs {'KIC'}];  end;
if job.ICs.KICc == 1, ICs = [ICs {'KICc'}]; end;

% execute operation
[N,M] = size(MS.SPMs);
for i = 1:N
    fprintf('\n-> Subject %d (%d out of %d):\n',i,i,N);
    for j = 1:M
        fprintf('   - Model %s (%d out of %d) ... ',MS.GLMs{j},j,M);
        load(MS.SPMs{i,j});                 % load SPM.mat
        MA_classic_ICs(SPM, [], ICs)        % calculate ICs
        fprintf('successful!\n');
    end;
end;
fprintf('\n');

% set output files
out = [];