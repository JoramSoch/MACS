function module = batch_MA_cvLME_auto
% _
% Configure MATLAB Batch for "MA: calculate cvLME (automatic)"
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 04:05 (V0.99/V15)
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

% Accuracy & Complexity
%-------------------------------------------------------------------------%
AnC        = cfg_menu;
AnC.tag    = 'AnC';
AnC.name   = 'Accuracy & Complexity';
AnC.val    = {0};
AnC.help   = {'Choose whether model accuracy and model complexity are calculated in addition to the log model evidence.'};
AnC.labels = {'No', 'Yes'};
AnC.values = {0, 1};


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MA: cvLME (auto)
%-------------------------------------------------------------------------%
module      = cfg_exbranch;
module.tag  = 'MA_cvLME_auto';
module.name = 'MA: calculate cvLME (automatic)';
module.val  = {MS_mat AnC};
module.help = {'Cross-Validated Log Model Evidence for General Linear Models'
               'Type "help MA_cvLME_multi" or "help MA_cvLME_single" for help.'};
module.prog = @run_module;


%=========================================================================%
% F U N C T I O N S                                                       %
%=========================================================================%

% Run batch
%-------------------------------------------------------------------------%
function out = run_module(job)

% get input variables
load(job.MS_mat{1});
AnC = logical(job.AnC);

% execute operation
[N,M] = size(MS.SPMs);
for i = 1:N
    fprintf('\n-> Subject %d (%d out of %d):\n',i,i,N);
    for j = 1:M
        fprintf('   - Model %s (%d out of %d) ... ',MS.GLMs{j},j,M);
        load(MS.SPMs{i,j});                 % load SPM.mat
        MA_cvLME_multi(SPM,[],[],AnC);      % calculate cvLME
        fprintf('successful!\n');
    end;
end;
fprintf('\n');

% set output files
out = [];