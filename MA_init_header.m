function H = MA_init_header(SPM, bin)
% _
% Initialize Header for Saving Brain Maps
% FORMAT H = MA_init_header(SPM, bin)
% 
%     SPM - a structure specifying an estimated GLM
%     bin - a logical indicating whether map is binary
% 
%     H   - a structure with header information
% 
% FORMAT H = MA_init_header(SPM, bin) generates a header structure from the
% information in SPM which is suited for binary data if bin is true and for
% real data if bin is false.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 24/10/2014, 18:10 (V0.2/V6)
%  Last edit: 18/08/2015, 07:40 (V0.3/V11)


% Initialize header
%-------------------------------------------------------------------------%
H = struct('fname',   '',...
           'dim',     SPM.VM.dim,...
           'dt',      [spm_type('float32') spm_platform('bigend')],...
           'mat',     SPM.VM.mat,...
           'pinfo',   SPM.VM.pinfo,...
           'descrip', '');

% Adapt for binary data
%-------------------------------------------------------------------------%
if bin
    H.dt = SPM.VM.dt;
end;