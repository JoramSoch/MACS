function MS_BMA_group(job, params, method)
% _
% Bayesian Model Averaging for General Linear Models (subject group)
% FORMAT MS_BMA_group(job, params, method)
%     job    - SPM batch editor "MS: perform BMA (manually)"
%     params - an M x P matrix of model parameters to be averaged
%     method - a string indicating analysis options (see below)
% 
% FORMAT MS_BMA_group(job, params, method) performs Bayesian model
% averaging according to the batch editor job with parameters of
% interest indexed by params and using the method indicated by method.
% 
% The input variable "params" can be an M x P matrix indexing parameters of
% interest for all models separately or a 1 x P vector if parameter indices
% are the same for all models.
% 
% The input variable "method" is a string containing up to five letters:
%     If it contains 'b', then the best model's parameters are extracted.
%     If it contains 'm', then the median model's parameters are extracted.
%     If it contains 'w', then the worst model's parameters are extracted.
%     If it contains 'r', then a random model's parameters are extracted.
%     If it contains 'a', then the averaged model parameters are estimated.
% 
% The default value of this input variable is 'ba' which means that the
% best model's parameters are extracted and BMA estimates are calculated.
% 
% Further information:
%     help MS_BMA_subject
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 03/03/2016, 18:35 (V0.4/V13)
%  Last edit: 11/05/2017, 18:20 (V1.0/V16)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;

% Get matlabbatch if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    design_mat = spm_select(1,'^*.mat','Select Batch Editor Job!');
    load(design_mat);
    job = matlabbatch{1}.spm.tools.MACS.MS_BMA_group_man;
    MS_BMA_group(matlabbatch);
    return
end;

% Set parameter indices if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(params), params = spm_input('parameter indices:',1,'r','[]'); end;

% Set analysis method if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(method), method = 'ba'; end;


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Get model parameters
%-------------------------------------------------------------------------%
N = numel(job.models);          % number of subjects
M = numel(job.models{1});       % number of models
S = 1; % already eliminated     % number of sessions
P = size(params,2);             % number of parameters

% Get log model evidence maps
%-------------------------------------------------------------------------%
LMEs = cell(N,M);
for i = 1:N
    for j = 1:M
        LMEs{i,j} = job.models{i}{j}{S};
    end;
end;

% Get model space directories
%-------------------------------------------------------------------------%
BMA_dirs = cell(N,1);
for i = 1:N
    BMA_dirs{i} = strcat(fileparts(LMEs{i,1}),'/','..');
end;

% Perform Bayesian model averaging
%-------------------------------------------------------------------------%
for i = 1:N
    fprintf('\n-> Subject %d (%d out of %d) ... ',i,i,N);
    MS_BMA_subject(LMEs(i,:)', params, method, BMA_dirs{i});
    fprintf('successful!\n');
end;
fprintf('\n');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);