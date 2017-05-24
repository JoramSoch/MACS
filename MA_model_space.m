function MA_model_space(job)
% _
% Specify Model Space of General Linear Models
% FORMAT MA_model_space(job)
%     job - SPM batch editor "MA: define model space"
% 
% FORMAT MA_model_space(job) specifies a model space consisting of M models
% applied to N subjects and saves all information into an MS.mat file.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 02:05 (V0.99/V15)
%  Last edit: 17/03/2017, 02:20 (V0.99/V15)


% Get model space dimensions
%-------------------------------------------------------------------------%
N = numel(job.models);          % number of subjects
M = numel(job.names);           % number of models

% Get working directory
%-------------------------------------------------------------------------%
MS.swd = job.dir{1};

% Assemble SPM.mats
%-------------------------------------------------------------------------%
MS.SPMs = cell(N,M);
for i = 1:N
    for j = 1:M
        MS.SPMs{i,j} = job.models{i}{j}{1};
    end;
end;

% Assemble GLM names
%-------------------------------------------------------------------------%
MS.GLMs = cell(1,M);
for j = 1:M
    MS.GLMs{j} = job.names{j};
end;

% Save model space file
%-------------------------------------------------------------------------%
save(strcat(MS.swd,'/','MS.mat'),'MS');