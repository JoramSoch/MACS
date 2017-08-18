% MACS Toolbox: Model Space Pipeline
% _
% Model Space Pipeline for the MACS Toolbox
% For details, see section 3.2 of the paper.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 18/08/2017, 17:35


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