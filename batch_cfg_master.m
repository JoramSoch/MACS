function module = batch_cfg_master
% _
% Configure MATLAB Batch Menu for MACS Toolbox (see spm_cfg.m)
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 17/03/2017, 00:45 (V0.99/V15)
%  Last edit: 09/02/2022, 11:40 (V1.4/V20)


% Define toolbox path
%-------------------------------------------------------------------------%
if ~isdeployed
    MACS_dir = dir(fullfile(spm('dir'),'toolbox','MACS*'));
    MACS_dir = MACS_dir([MACS_dir.isdir]);
    if ~strcmp(MACS_dir(1).name,'MACS')
        old_dir = fullfile(spm('dir'),'toolbox',MACS_dir(1).name);
        new_dir = fullfile(spm('dir'),'toolbox','MACS');
        movefile(old_dir, new_dir, 'f');
    end;
    addpath(fullfile(spm('dir'),'toolbox','MACS'));
end;


%=========================================================================%
% M O D U L E                                                             %
%=========================================================================%

% MACS Toolbox
%-------------------------------------------------------------------------%
module        = cfg_choice;
module.tag    = 'MACS';
module.name   = 'MACS Toolbox';
module.help   = {'Model assessment, comparison and selection (MACS)'};
module.values = {batch_MA_model_space, ...      % MA: define model space
                 batch_MA_inspect_GoF, ...      % MA: inspect goodness of fit
                 batch_MA_classic_ICs_man, ...  % MA: classical ICs (manually)
                 batch_MA_classic_ICs_auto, ... % MA: classical ICs (automatic)
                 batch_MA_cvLME_man, ...        % MA: calculate cvLME (manually)
                 batch_MA_cvLME_auto, ...       % MA: calculate cvLME (automatic)
                 batch_MA_LFE_group_man, ...    % MA: calculate LFE (manually)
                 batch_MA_LFE_group_auto, ...   % MA: calculate LFE (automatic)
                 batch_MC_LBF_group_man, ...    % MC: calculate LBF (manually)
                 batch_MC_LBF_group_auto, ...   % MC: calculate LBF (automatic)
                 batch_MS_PPs_group_man, ...    % MS: calculate PPs (manually)
                 batch_MS_PPs_group_auto, ...   % MS: calculate PPs (automatic)
                 batch_MS_BMS_group_man, ...    % MS: perform BMS (manually)
                 batch_MS_BMS_group_auto, ...   % MS: perform BMS (automatic)
                 batch_MS_BMA_group_man, ...    % MS: perform BMA (manually)
                 batch_MS_BMA_group_auto, ...   % MS: perform BMA (automatic)
                 batch_MS_BMS_fams_man, ...     % MS: perform family BMS (manually)
                 batch_MS_BMS_fams_auto, ...    % MS: perform family BMS (automatic)
                 batch_MS_SMM_BMS, ...          % MS: generate SMM from BMS
                 batch_MF_visualize             % MF: visualize high-dimensional data
                };