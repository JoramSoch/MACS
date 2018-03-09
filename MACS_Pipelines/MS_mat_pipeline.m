% MACS Toolbox: Model Space Pipeline
% _
% This script assists in setting up a batch for defining a model space that
% can be viewed and executed in the SPM batch editor. It is particularly
% advantageous when the number of subjects or the number of models in your
% analyses is very large.
% 
% This script assumes that you have organized your data in the form of a
% subject-model hierarchy looking like this:
% 
%     [stat_dir]\
%     [work_dir]\
%         sub01\
%             mod01\
%             mod02\
%             mod03\
%             ...
%             mod08\
%             mod09\
%             mod10\
%         sub02\
%             mod01\
%             ...
%             mod10\
%         sub03\
%         ...
%         sub23\
%         sub24\
%         sub25\
% 
% If this is the case, you can simply enter
% - the statistics directory into "stat_dir",
% - the working directory into "work_dir",
% - the subject folder names into "subj_ids" and
% - the model folder names into "mod_nams" below.
% 
% In addition, you will have to specify
% - a model space name as "ms_name" and
% - a model space suffix as "ms_suff"
% which, together with the statistics directory, will determine where the
% model space directory will be located and the model space file will be
% written. Use these two parameters to distinguish different model space
% and analyses from each other.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/08/2017, 17:35 (V1.1/V17)
%  Last edit: 09/03/2018, 09:30 (V1.2/V18)


%%% Step 0: Study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project directories
stat_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\analyses';
work_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\functional';
% Note: This was the data (Bogler et al., 2013) analyzed in the paper (Soch et al., 2016).

% list of subjects
subj_ids = {'sub01' 'sub02' 'sub03' 'sub04' 'sub05' 'sub06' 'sub07' 'sub08' 'sub09' 'sub10' ...
            'sub11' 'sub12' 'sub13' 'sub14' 'sub15' 'sub16' 'sub17' 'sub18' 'sub19' 'sub20' ...
            'sub21' 'sub22' 'sub23' 'sub24' 'sub25'};

% list of models
mod_nams = {'mod01' 'mod02' 'mod03' 'mod04' 'mod05' 'mod06' 'mod07' 'mod08' 'mod09' 'mod10'};

% model space details
ms_name  =  'MS01';
ms_suff  =  'test';

% study dimensions
N = numel(subj_ids);
M = numel(mod_nams);


%%% Step 1: Create model space job %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% working directory
job.dir{1} = strcat(stat_dir,'/',ms_name,'_',ms_suff,'/');

% assemble SPM.mats
for i = 1:N
    for j = 1:M
        job.models{i}{j}{1} = strcat(work_dir,'/',subj_ids{i},'/',mod_nams{j},'/','SPM.mat');
    end;
end;

% assemble GLM names
for j = 1:M
    job.names{j} = mod_nams{j};
end;


%%% Step 2: Execute model space job %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save batch
clear matlabbatch
mkdir(job.dir{1});
filename = strcat(job.dir{1},'batch.mat');
matlabbatch{1}.spm.tools.MACS.MA_model_space = job;
save(filename,'matlabbatch');

% execute job
MA_model_space(job);

% display message
fprintf('\n');
fprintf('\n-> Thank you! The following files have been created:\n');
fprintf('   - SPM batch: %s.\n', strcat(job.dir{1},'batch.mat'));
fprintf('   - MS.mat file: %s.\n', strcat(job.dir{1},'MS.mat'));
fprintf('\n');