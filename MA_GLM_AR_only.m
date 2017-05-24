function MA_GLM_AR_only(SPM_mat)
% _
% General Linear Model with AR Component Estimation only
% FORMAT MA_GLM_AR_only(SPM_mat)
%     SPM_mat - a string indicating a specified GLM
% 
% FORMAT MA_GLM_AR_only(SPM_mat) estimates a general linear model and
% saves only the auto-regressive components used for non-sphericity
% correction. These are estimated from the residual auto-correlations using
% a restricted maximum likelihood (ReML) approach and saved into SPM.xVi.V.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 18/08/2015, 04:10 (V0.3/V11)
%  Last edit: 18/08/2015, 08:30 (V0.3/V11)


% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
end;

% Load SPM.mat if available
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    load(SPM_mat);
    SPM.swd = fileparts(SPM_mat);
    cd(SPM.swd);
catch
    MA_GLM_AR_only;
    return
end;

% If current version is SPM8
%-------------------------------------------------------------------------%
if strcmp(spm('Ver'),'SPM8')
    
    % Estimate general linear model
    %---------------------------------------------------------------------%
    SPM = spm_spm(SPM);
    
    % Delete beta images
    %---------------------------------------------------------------------%
    for i = 1:length(SPM.Vbeta)
        filename = strcat(SPM.swd,'/',SPM.Vbeta(i).fname);
        [SPM_dir, name, ext] = fileparts(filename);
        if strcmp(ext,'.img')       % SPM8 or earlier
            spm_unlink(strcat(SPM.swd,'/',name,'.hdr'));
            spm_unlink(strcat(SPM.swd,'/',name,'.img'));
        end;
        if strcmp(ext,'.nii')       % SPM12 or later
            spm_unlink(strcat(SPM.swd,'/',name,'.nii'));
        end;
    end;
    
    % Delete contrast maps
    %---------------------------------------------------------------------%
    if ~isempty(SPM.xCon)
        for i = 1:length(SPM.xCon)
            for j = 1:2
                if j == 1, filename = strcat(SPM.swd,'/',SPM.xCon(i).Vcon.fname); end;
                if j == 2, filename = strcat(SPM.swd,'/',SPM.xCon(i).Vspm.fname); end;
                [SPM_dir, name, ext] = fileparts(filename);
                if strcmp(ext,'.img')       % SPM8 or earlier
                    spm_unlink(strcat(SPM.swd,'/',name,'.hdr'));
                    spm_unlink(strcat(SPM.swd,'/',name,'.img'));
                end;
                if strcmp(ext,'.nii')       % SPM12 or later
                    spm_unlink(strcat(SPM.swd,'/',name,'.nii'));
                end;
            end;
        end;
    end;
    
end;

% If current version is SPM12
%-------------------------------------------------------------------------%
if strcmp(spm('Ver'),'SPM12')
    
    % Estimate non-sphericity
    %---------------------------------------------------------------------%
    [xVi, M] = spm_est_non_sphericity(SPM);
    
    % Insert non-sphericity
    %---------------------------------------------------------------------%
    W = spm_sqrtm(spm_inv(xVi.V));
    W = W.*(abs(W) > 1e-6);
    SPM.xX.W = W;
    SPM.xVi  = xVi;
    
    % Insert masking image
    %---------------------------------------------------------------------%
    H = SPM.xY.VY(1);
    H.fname   = 'mask.nii';
    H.dt      = [spm_type('uint8') spm_platform('bigend')];
    H.pinfo   = [1 0 0]';
    H.descrip = 'MA_GLM_AR_only: masking image';
    H = spm_data_hdr_write(H);
    spm_write_vol(H,double(M));
    SPM.VM = H;
    
    % Save SPM.mat
    %---------------------------------------------------------------------%
    save(SPM_mat,'SPM');
    
end;

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);