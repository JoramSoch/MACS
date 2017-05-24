function MS_BMA_subject(LMEs, para, opts, folder)
% _
% Bayesian Model Averaging for General Linear Models (single subject)
% FORMAT MS_BMA_subject(LMEs, para, opts, folder)
%     LMEs   - an M x 1 cell array specifying log model evidence maps
%     para   - an M x P matrix of model parameters to be averaged
%     opts   - a string indicating analysis options (see below)
%     folder - a string indicating the folder for the output images
% 
% FORMAT MS_BMA_subject(LMEs, para, opts, folder) performs Bayesian model
% averaging based on the log model evidence maps LMEs using the settings
% indicated by opts and saves averaged parameters indexed by para into
% (sub-directories of) the directory folder.
% 
% The model evidence is the probability of the data y (e.g. fMRI signal),
% just given the model m (e.g. a GLM), regardless of any particular
% parameter values b (e.g. regression weights) ([1], eq. 3):
%     p(y|m) = INT p(y|b,m) p(b|m) db
% Then, the model evidence can be used to calculate posterior probabilties
% (PP) using Bayes' Theorem ([1], eq. 2):
%     p(m_i|y) = [p(y|m_i) * p(m_i)] / SUM_j^M [p(y|m_j) p(m_j)]
% Then, these PPs can be used to calculate averaged model parameters as
% weighted averages ([1], eq. 1):
%     b_BMA = SUM_i^M [b_i * p(m_i|y)]
% These are the three steps of Bayesian model averaging (BMA). In fMRI,
% BMA has been implemented for dynamic causal models (DCMs) [2] and now
% also for general linear models (GLMs) [3]. In our implementation, we
% use the cross-validated log model evidence (cvLME) to calculate PPs.
% 
% The input variable "LMEs" is an M x 1 cell array indicating filenames of
% cvLME maps, where M is the number of models.
% 
% The input variable "para" can be an M x P matrix indexing parameters of
% interest for all models separately or a 1 x P vector if parameter indices
% are the same for all models, where P is the number of parameters.
% 
% The input variable "opts" is a string which contains up to five letters:
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
%     help MS_BMA_group
% 
% References:
% [1] Hoeting JA, Madigan D, Raftery AE, Volinsky CT (1999):
%     "Bayesian Model Averaging: A Tutorial".
%     Statistical Science, vol. 14, no. 4, pp. 382-417.
% [2] Penny WD, Stephan KE, Daunizeau J, Rosa MJ, Friston KJ, Schofield TM,
%     Leff AP (2010): "Comparing Families of Dynamic Causal Models".
%     PLoS ONE, vol. 6, iss. 3, e1000709.
% [3] Soch J, Meyer AP, Haynes JD, Allefeld C (2017): "How to improve parameter estimates in 
%     GLM-based fMRI data analysis: cross-validated Bayesian model averaging".
%     NeuroImage, in review. URL: http://biorxiv.org/content/early/2016/12/20/095778
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 03/03/2016, 13:30 (V0.4/V13)
%  Last edit: 12/05/2017, 13:25 (V1.0/V16)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get current directory
%-------------------------------------------------------------------------%
orig_dir = pwd;

% Get model parameters
%-------------------------------------------------------------------------%
M = size(LMEs,1);               % number of models
P = size(para,2);               % number of parameters

% Get image dimensions
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});           % LME image header
V = prod(H.dim);                % number of voxels

% Expand parameters if necessary
%-------------------------------------------------------------------------%
if size(para,1) == 1, para = repmat(para,[M 1]); end;

% Set options if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(opts), opts = 'ba'; end;

% Set analysis directory
%-------------------------------------------------------------------------%
BMA_dir = strcat(folder,'/','MS_BMA_subject','_',opts);

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: load');

% Load log model evidences
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Load log model evidences...' , '');
LME = zeros(M,V);               % M x V matrix of LMEs
for i = 1:M                     % models
    lme_hdr = spm_vol(LMEs{i});
    lme_img = spm_read_vols(lme_hdr);
    lme_img = reshape(lme_img,[1 V]);
    LME(i,:)= lme_img;
    % NOTE: This subject-level BMA operates on model-wise, but
    % not session-wise LMEs, because session-wise (oos)LMEs have
    % already been summed up during first-level model assessment.    
    spm_progress_bar('Set',(i/M)*100);
end;
clear lme_hdr lme_img

% Load parameter estimates
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Load parameter estimates...' , '');
B = zeros(M,V,P);               % M x V x P array of betas
for i = 1:M                     % models
    SPM_dir = fileparts(LMEs{i});
    SPM_mat = strcat(SPM_dir,'/','SPM.mat');
    load(SPM_mat);
    cd(SPM.swd);
    S = numel(SPM.Sess);
    for j = 1:P                 % parameters
        for k = 1:S             % sessions
            if para(i,j) ~= 0   % indexed parameter
                beta_hdr = SPM.Vbeta(SPM.Sess(k).col(para(i,j)));
                beta_img = spm_read_vols(beta_hdr);
                beta_img = reshape(beta_img,[1 V]);
                B(i,:,j) = B(i,:,j) + beta_img;
            else                % zero parameter
                B(i,:,j) = B(i,:,j) + zeros(1,V);
            end;
            spm_progress_bar('Set',(((i-1)*P*S+(j-1)*S+k)/(M*P*S))*100);
        end;                    % average
        B(i,:,j) = 1/S * B(i,:,j);
    end;
end;
clear SPM_dir SPM_mat beta_hdr beta_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Create mask image
%-------------------------------------------------------------------------%
[m_img m_hdr m_ind] = MS_create_mask(LME, H);


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: estimate');

% Select in-mask voxels only
%-------------------------------------------------------------------------%
v = length(m_ind);
d = floor(v/100);

% Get best model's parameters
%-------------------------------------------------------------------------%
if ismember('b',opts)
    spm_progress_bar('Init', 100, 'Get best model''s parameters...' , '');
    Bb    = NaN(P,V);           % best model's betas
    [c,i] = max(LME,[],1);
    for j = 1:P                 % parameters
        for l = 1:v             % voxels
            Bb(j,m_ind(l)) = B(i(m_ind(l)),m_ind(l),j);
        end;
        spm_progress_bar('Set',(j/P)*100);
    end;
end;

% Get median model's parameters
%-------------------------------------------------------------------------%
if ismember('m',opts)
    spm_progress_bar('Init', 100, 'Get median model''s parameters...' , '');
    Bm    = NaN(P,V);           % median model's betas
    for j = 1:P                 % parameters
        for l = 1:v             % voxels
            [c,i] = min(abs(LME(:,m_ind(l),j)-median(LME(:,m_ind(l),j))));
            Bm(j,m_ind(l)) = B(i(1),m_ind(l),j);
        end;
        spm_progress_bar('Set',(j/P)*100);
    end;
end;

% Get worst model's parameters
%-------------------------------------------------------------------------%
if ismember('w',opts)
    spm_progress_bar('Init', 100, 'Get worst model''s parameters...' , '');
    Bw    = NaN(P,V);           % worst model's betas
    [c,i] = min(LME,[],1);
    for j = 1:P                 % parameters
        for l = 1:v             % voxels
            Bw(j,m_ind(l)) = B(i(m_ind(l)),m_ind(l),j);
        end;
        spm_progress_bar('Set',(j/P)*100);
    end;
end;

% Get random model's parameters
%-------------------------------------------------------------------------%
if ismember('r',opts)
    spm_progress_bar('Init', 100, 'Get random model''s parameters...' , '');
    Br    = NaN(P,V);           % random model's betas
    i     = randi([1 M],[1 V]);
    for j = 1:P                 % parameters
        for l = 1:v             % voxels
            Br(j,m_ind(l)) = B(i(m_ind(l)),m_ind(l),j);
        end;
        spm_progress_bar('Set',(j/P)*100);
    end;    
end;

% Get averaged model parameters
%-------------------------------------------------------------------------%
if ismember('a',opts)
    spm_progress_bar('Init', 100, 'Get averaged model parameters...' , '');
    Ba    = NaN(P,V);           % averaged betas
    % obtain posterior model probabilities
    prior = 1/M * ones(M,1);
    LMEp  = LME - repmat(mean(LME,1),[M 1]);
    LMEp  = exp(LMEp) .* repmat(prior,[1 V]);
    post  = LMEp ./ repmat(sum(LMEp,1),[M 1]);
    % obtain averaged parameter estimates
    for j = 1:P                 % parameters
        Ba(j,:) = sum(B(:,:,j).*post,1);
        spm_progress_bar('Set',(j/P)*100);
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','MS_BMA_subject: save');

% Initialise image files
%-------------------------------------------------------------------------%
H = spm_vol(LMEs{1});
if ~exist(BMA_dir,'dir')
    mkdir(BMA_dir);
end;

% Save best model's parameters
%-------------------------------------------------------------------------%
if ismember('b',opts)
    cd(BMA_dir);
    mkdir('best_model');
    cd('best_model');
    for j = 1:P                 % parameters
        H.fname   = strcat('beta_',MF_int2str0(para(1,j),4),'.nii');
        H.descrip = 'MS_BMA_subject: best model''s parameters';
        spm_write_vol(H,reshape(Bb(j,:),H.dim));
    end;
end;

% Save median model's parameters
%-------------------------------------------------------------------------%
if ismember('m',opts)
    cd(BMA_dir);
    mkdir('median_model');
    cd('median_model');
    for j = 1:P                 % parameters
        H.fname   = strcat('beta_',MF_int2str0(para(1,j),4),'.nii');
        H.descrip = 'MS_BMA_subject: median model''s parameters';
        spm_write_vol(H,reshape(Bm(j,:),H.dim));
    end;
end;

% Save worst model's parameters
%-------------------------------------------------------------------------%
if ismember('w',opts)
    cd(BMA_dir);
    mkdir('worst_model');
    cd('worst_model');
    for j = 1:P                 % parameters
        H.fname   = strcat('beta_',MF_int2str0(para(1,j),4),'.nii');
        H.descrip = 'MS_BMA_subject: worst model''s parameters';
        spm_write_vol(H,reshape(Bw(j,:),H.dim));
    end;
end;

% Save random model's parameters
%-------------------------------------------------------------------------%
if ismember('r',opts)
    cd(BMA_dir);
    mkdir('random_model');
    cd('random_model');
    for j = 1:P                 % parameters
        H.fname   = strcat('beta_',MF_int2str0(para(1,j),4),'.nii');
        H.descrip = 'MS_BMA_subject: random model''s parameters';
        spm_write_vol(H,reshape(Br(j,:),H.dim));
    end;
end;

% Save averaged model parameters
%-------------------------------------------------------------------------%
if ismember('a',opts)
    cd(BMA_dir);
    for j = 1:P                 % parameters
        H.fname   = strcat('beta_',MF_int2str0(para(1,j),4),'_BMA','.nii');
        H.descrip = 'MS_BMA_subject: averaged model parameters';
        spm_write_vol(H,reshape(Ba(j,:),H.dim));
    end;
end;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);