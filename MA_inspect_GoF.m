function varargout = MA_inspect_GoF(varargin)
% _
% Inspect Goodness of Fit in a General Linear Model
% FORMAT varargout = MA_inspect_GoF(varargin)
% 
%     varargin  - cell array of input variables
%     varargout - cell array of output variables
% 
% FORMAT SPM = MA_inspect_GoF('Specify') opens SPM's file selector and lets
% you select an SPM.mat representing an estimated general linear model.
% 
%     SPM - structure specifying an estimated GLM; important fields:
%     o swd   - statistics working directory
%     o VM    - header information for mask file
%     o VY    - header information for data files
%     o xX    - information on the design matrix
%     o Vbeta - header information for beta image
%     o Sess  - information on the different runs
% 
% FORMAT GLM = MA_inspect_GoF('LoadData',SPM) uses an SPM structure and
% loads data, design and parameter estimates for that estimated model.
% 
%     GLM - structure specifying an estimated GLM; important fields:
%     o Y       - n x v data matrix
%     o X       - n x p design matrix
%     o B       - p x v parameter matrix
%     o K, W, V - filtering, whitening, non-sphericity matrix
%     o n, p, v - number of data points, parameters, in-mask voxels
%     o w       - total number of voxels including out-mask voxels
%     o M       - 1 x w binary mask image
%     o XYZ     - 3 x w coordinate matrix
%     o KWX     - whitened and filtered design
%     o KWY     - whitened and filtered signal
%     o KWZ     - whitened and filtered prediction
%     o KWXc    - baseline-corrected design
%     o KWYc    - baseline-corrected signal
%     o KWZc    - baseline-corrected prediction
%     o pc      - non-baseline regressor indices
%     o sig2    - 1 x v residual variances
%     o R2      - 1 x v coefficients of determination
%     o adj_R2  - 1 x v adjusted R^2 values
%     o gen_R2  - 1 x v generalized R^2 values
%     o mf_SNR  - 1 x v model-free SNR values
%     o mb_SNR  - 1 x v model-based SNR values
% 
% FORMAT [hReg,xSPM] = MA_inspect_GoF('Setup',SPM) takes an SPM.mat, loads
% the corresponding data and prepares figures, graphics and handles for
% data visualization. The function should be called in this mode when SPM
% is loaded externally and not using SPM's file selector.
% 
% For several reasons, the implicit baseline is subtracted from the
% (whitened and filtered) measured as well as from the predicted signal.
% First of all, the absolute value of the baseline is meaningless and
% arbitrary. Second, the baseline differs across sessions which can make it
% hard to inspect the signal when between-session differences are large.
% Third, there are heavy spikes at the start and end of each session when
% the signal is whitened which can be avoided when the (whitened) baseline
% is subtracted.
% 
% This function was written by modifying the SPM Results User Interface.
% Consequently, much of this code is stolen from spm_results_ui.m!
% Additional inspiration came from spm_run_bms_vis.m.
% 
% References:
% [1] Soch J, Allefeld C (2018): "MACS - a new SPM toolbox for model
%     assessment, comparison and selection". Journal of Neuroscience
%     Methods, vol. 306, pp. 19-31.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 12/02/2015, 05:35 (V0.3/V10)
%  Last edit: 09/03/2018, 11:40 (V1.2/V18)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Declare persistent variables
%-------------------------------------------------------------------------%
persistent SPM GLM orig_dir
% SPM      - structure describing an estiamted GLM
% GLM      - structure with data from an estimated GLM
% orig_dir - working directory when setup mode is called

% Set action if required
%-------------------------------------------------------------------------%
if nargin == 0
    action = 'Specify';
else
    action = varargin{1};
end

% Check action string
%-------------------------------------------------------------------------%
switch lower(action)


%=========================================================================%
% Action: 'Specify'                                 specify visualization %
%=========================================================================%
case 'specify'
% SPM = MA_inspect_GoF('Specify')

    % Initialise
    %---------------------------------------------------------------------%
    Finter = spm('FigName','MA_inspect_GoF: Specify');
    
    % Select and load SPM.mat
    %---------------------------------------------------------------------%
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    load(SPM_mat);
    
    % Select type of plot
    %---------------------------------------------------------------------%
    PlotTypes    = {'linegraph','scatterplot'};
    SPM.PlotType = PlotTypes{spm_input('Plot Type:',1,'b',PlotTypes,[1 2])};

    % Set up visualization
    %---------------------------------------------------------------------%
    MA_inspect_GoF('Setup',SPM);

    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {SPM};
    

%=========================================================================%
% Action: 'Setup'                                    set up visualization %
%=========================================================================%
case 'setup'
% [hReg,xSPM] = MA_inspect_GoF('Setup',SPM)

    % Return
    %---------------------------------------------------------------------%
    if nargin < 2
        SPM = MA_inspect_GoF('Specify');
        return
    end;
    
    % Initialise
    %---------------------------------------------------------------------%
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MA_inspect_GoF: Setup');
    spm_clf('Satellite');
    
    % Change directory
    %---------------------------------------------------------------------%
    orig_dir = pwd;
 
    % Load input images
    %---------------------------------------------------------------------%
    SPM = varargin{2}; cd(SPM.swd);
    GLM = MA_inspect_GoF('LoadData',SPM);

    % Expand to all voxels
    %---------------------------------------------------------------------%    
    R2_all  = NaN(1,GLM.w);
    R2_all(GLM.m) = GLM.R2;

    % Find activated voxels
    %---------------------------------------------------------------------%
    thresh  = MF_quantile(R2_all(~isnan(R2_all)),[0.05 0.95]);
    stv_ind = find(R2_all < thresh(1) | R2_all > thresh(2));
    XYZ_ord      = zeros(3,GLM.w);          % ordinal voxel co-ordinates
    XYZ_ord(1,:) = kron(ones(1,SPM.VM.dim(3)*SPM.VM.dim(2)), [1:SPM.VM.dim(1)]);
    XYZ_ord(2,:) = kron(ones(1,SPM.VM.dim(3)), kron([1:SPM.VM.dim(2)], ones(1,SPM.VM.dim(1))));
    XYZ_ord(3,:) = kron([1:SPM.VM.dim(3)], ones(1,SPM.VM.dim(2)*SPM.VM.dim(1)));
    XYZ_ind      = XYZ_ord(:,stv_ind);      % supra-threshold voxel indexes
    
    % Create xSPM structure
    %---------------------------------------------------------------------%
    xSPM.swd       = SPM.swd;
    xSPM.title     = 'Goodness of Fit';
    xSPM.Z         = R2_all(stv_ind);
    xSPM.M         = SPM.VM.mat;
    xSPM.iM        = inv(xSPM.M);
    xSPM.VOX       = sqrt(sum(xSPM.M(1:3,1:3).^2));
    xSPM.DIM       = SPM.VM.dim(1:3)'; 
    xSPM.XYZ       = XYZ_ind;
    xSPM.XYZmm     = xSPM.M * [xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    xSPM.XYZmm     = xSPM.XYZmm(1:3,:);
    xSPM.units     = {'mm' 'mm' 'mm'};
    xSPM.thresDesc = 'none';
    
    % Get space information
    %---------------------------------------------------------------------%
    M   = xSPM.M;
    DIM = xSPM.DIM;
    
    % Set up visualize GUI
    %---------------------------------------------------------------------%
    spm_clf(Finter);
    spm('FigName',['MA_inspect_GoF: Setup'],Finter);
    hReg = MA_inspect_GoF('SetupGUI',M,DIM,xSPM,Finter);

    % Activate section overlay
    %---------------------------------------------------------------------%
    section = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
    spm_sections(xSPM,hReg,section);
 
    % Store handles of graphics window objects
    %---------------------------------------------------------------------%
    h = get(Fgraph,'Children');
    h = findobj(h,'flat','HandleVisibility','on');
    h = findobj(h);
    g = get(h,'Visible');
    % set(hResAx,'Tag','PermRes','UserData',struct('h',h,'g',{g}))
 
    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {hReg,xSPM};
    spm('Pointer','Arrow');


%=========================================================================%
% Action: 'LoadData'                                    load input images %
%=========================================================================%
case 'loaddata'
% GLM = MA_inspect_GoF('LoadData',SPM)
    
    % Get model parameters
    %---------------------------------------------------------------------%
    GLM.X = SPM.xX.X;           % design matrix
    GLM.K = SPM.xX.K;           % filtering matrix
    GLM.W = SPM.xX.W;           % whitening matrix
    GLM.V = SPM.xVi.V;          % non-sphericity
    GLM.n = size(GLM.X,1);      % number of data points
    GLM.p = size(GLM.X,2);      % number of parameters
    GLM.d = ceil(GLM.n/100);
    
    % Load mask image
    %---------------------------------------------------------------------%
    [GLM.M m_dim GLM.m] = MA_load_mask(SPM);
    [m_img GLM.XYZ] = spm_read_vols(SPM.VM);
    GLM.w = prod(SPM.VM.dim);   % number of all voxels
    GLM.v = length(GLM.m);      % number of in-mask voxels
    clear m_img m_dim
    
    % Load time series
    %---------------------------------------------------------------------%
    spm_progress_bar('Init',100,'Load time series...','');
    GLM.Y = zeros(GLM.n,GLM.v);
    for i = 1:GLM.n
        y_img = spm_read_vols(SPM.xY.VY(i));
        y_img = reshape(y_img,[1 GLM.w]);
        GLM.Y(i,:) = y_img(GLM.m);
        if mod(i,GLM.d) == 0, spm_progress_bar('Set',(i/GLM.n)*100); end;
    end;
    clear y_img
    spm_progress_bar('Clear');
    
    % Load parameter estimates
    %---------------------------------------------------------------------%
    spm_progress_bar('Init',100,'Load parameter estimates...','');
    GLM.B = zeros(GLM.p,GLM.v);
    for j = 1:GLM.p
        b_img = spm_read_vols(SPM.Vbeta(j));
        b_img = reshape(b_img,[1 GLM.w]);
        GLM.B(j,:) = b_img(GLM.m);
        spm_progress_bar('Set',(j/GLM.p)*100);
    end;
    clear b_img
    spm_progress_bar('Clear');
    
    % Measured and predicted signal
    %---------------------------------------------------------------------%
    GLM.KWX  = spm_filter(GLM.K,GLM.W*GLM.X);
    GLM.KWY  = spm_filter(GLM.K,GLM.W*GLM.Y);
    GLM.KWZ  = GLM.KWX * GLM.B;
    GLM      = rmfield(GLM,{'Y'});
    
    % Correct for implicit baseline
    %---------------------------------------------------------------------%
    GLM.pc   = [1:(SPM.xX.iB(1)-1)];
    GLM.base = GLM.KWX(:,SPM.xX.iB) * GLM.B(SPM.xX.iB,:);
    GLM.KWXc = GLM.KWX(:,GLM.pc);
    GLM.KWYc = GLM.KWY - GLM.base;
    GLM.KWZc = GLM.KWZ - GLM.base;
    GLM      = rmfield(GLM,{'KWX', 'KWY', 'KWZ', 'base'});
    
    % Calculate GoF as overlay
    %---------------------------------------------------------------------%
    [GLM.sig2, GLM.R2, GLM.adj_R2, GLM.gen_R2] = ME_GLM_GoF(GLM.KWYc, GLM.KWXc, [], GLM.B(GLM.pc,:));
    [GLM.mf_SNR, GLM.mb_SNR] = ME_GLM_SNR(GLM.KWYc, GLM.KWXc, [], GLM.B(GLM.pc,:));
    
    % Save GoF as images
    %---------------------------------------------------------------------%
    H = MA_init_header(SPM, false);
    GoF_img = NaN(size(GLM.M));
    GoF_img(GLM.m) = GLM.R2;                % R^2
    H.fname   = 'MA_GoF_R2.nii';
    H.descrip = 'MA_inspect_GoF: coefficient of determination (R^2)';
    spm_write_vol(H,reshape(GoF_img,SPM.VM.dim));
    SPM.MACS.R2 = H;
    GoF_img(GLM.m) = GLM.adj_R2;            % adj. R^2
    H.fname   = 'MA_GoF_R2_adj.nii';
    H.descrip = 'MA_inspect_GoF: adjusted coefficient of determination (adj. R^2)';
    spm_write_vol(H,reshape(GoF_img,SPM.VM.dim));
    SPM.MACS.R2_adj = H;
    GoF_img(GLM.m) = GLM.mf_SNR;            % mf. SNR
    H.fname   = 'MA_GoF_SNR_mf.nii';
    H.descrip = 'MA_inspect_GoF: model-free signal-to-noise ratio (mf. SNR)';
    spm_write_vol(H,reshape(GoF_img,SPM.VM.dim));
    SPM.MACS.SNR_mf = H;
    GoF_img(GLM.m) = GLM.mb_SNR;            % mb. SNR
    H.fname   = 'MA_GoF_SNR_mb.nii';
    H.descrip = 'MA_inspect_GoF: model-based signal-to-noise ratio (mb. SNR)';
    spm_write_vol(H,reshape(GoF_img,SPM.VM.dim));
    SPM.MACS.SNR_mb = H;
    save(strcat(SPM.swd,'/','SPM.mat'),'SPM');
    clear H GoF_img
    
    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {GLM};


%=========================================================================%
% Action: 'SetupGUI'                                 set up visualize GUI %
%=========================================================================%
case 'setupgui'
% hReg = MA_inspect_GoF('SetupGUI',M,DIM,xSPM,Finter)

    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 5, Finter = 'Interactive'; else Finter = varargin{5}; end
    if nargin < 4, error('Insufficient arguments'); end
    M      = varargin{2};
    DIM    = varargin{3};
    Finter = spm_figure('GetWin',Finter);
    WS     = spm('WinScale');
    FS     = spm('FontSizes');

    % Create frame for GUI objects
    %---------------------------------------------------------------------%
    hPan = uipanel('Parent',Finter,'Title','','Units','Pixels',...
            'Position',[001 001 400 190].*WS,...
            'BorderType','Line', 'HighlightColor',[0 0 0],...
            'BackgroundColor',spm('Colour'));
    hReg = uipanel('Parent',hPan,'Title','','Units','Pixels',...
            'BorderType','Etchedin', ...
            'Position',[005 005 390 180].*WS,...
            'BackgroundColor',[179 179 179]/255);

    % Initialise registry in hReg frame object
    %---------------------------------------------------------------------%
    [hReg,xyz] = spm_XYZreg('InitReg',hReg,M,DIM,[0;0;0]);

    % Setup editable XYZ widgets & cross-register
    %---------------------------------------------------------------------%
    hFxyz      = MA_inspect_GoF('DrawXYZgui',M,DIM,varargin{4},xyz,hReg);
    spm_XYZreg('XReg',hReg,hFxyz,'MA_inspect_GoF');

    % Setup buttons for visualization functions
    %---------------------------------------------------------------------%
    MA_inspect_GoF('DrawButts',hReg,DIM,Finter,WS,FS);
    set(findobj(hPan),'Units','Normalized','FontUnits','Normalized');

    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {hReg};
 

%=========================================================================%
% Action: 'DrawButts'                  draw buttons in interactive window %
%=========================================================================%
case 'drawbutts'
% MA_inspect_GoF('DrawButts',hReg,DIM,Finter,WS,FS)

    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 3, error('Insufficient arguments'); end
    hReg = varargin{2};
    DIM  = varargin{3};
    if nargin < 4, Finter = spm_figure('FindWin','Interactive'); else Finter = varargin{4}; end
    if nargin < 5, WS = spm('WinScale');  else WS = varargin{5}; end
    if nargin < 6, FS = spm('FontSizes'); else FS = varargin{6}; end

    % Draw GUI control buttons
    %---------------------------------------------------------------------%
    uicontrol('Parent',hReg,'Style','PushButton','String','clear',...
        'ToolTipString','Clear results subpane',...
        'FontSize',FS(9),'ForegroundColor','b',...
        'Callback',['MA_inspect_GoF(''Clear''); ',...
          'spm_input(''!DeleteInputObj''),',...
          'spm_clf(''Satellite'')'],...
        'Interruptible','on','Enable','on',...
        'DeleteFcn','spm_clf(''Graphics'')',...
        'Position',[280 050 048 020].*WS);
    uicontrol('Parent',hReg,'Style','PushButton','String','exit',...
        'ToolTipString','Exit the results section',...
        'FontSize',FS(9),'ForegroundColor','r',...
        'Callback','MA_inspect_GoF(''close'')',...
        'Interruptible','on','Enable','on',...
        'Position',[332 050 048 020].*WS);

    
%=========================================================================%
% Action: 'DrawXYZgui'                 draw XYZ GUI in interactive window %
%=========================================================================%
case 'drawxyzgui'
% hFxyz = MA_inspect_GoF('DrawXYZgui',M,DIM,xSPM,xyz,hReg)

    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 6, hReg = spm_XYZreg('FindReg','Interactive'); else hReg = varargin{6}; end
    if nargin < 5, xyz = [0;0;0]; else xyz = varargin{5}; end
    if nargin < 4, error('Insufficient arguments'); end
    DIM = varargin{3};
    M   = varargin{2};
    xyz = spm_XYZreg('RoundCoords',xyz,M,DIM);
    WS  = spm('WinScale');
    FS  = spm('FontSizes');
    PF  = spm_platform('fonts');

    % Create XYZ control objects
    %---------------------------------------------------------------------%
    hFxyz = uipanel('Parent',hReg,'Title','co-ordinates','Units','Pixels',...
        'Position',[005 005 265 040].*WS,...
        'BorderType','Beveledout',...
        'ShadowColor',[0.5 0.5 0.5],...
        'FontAngle','Italic',...
        'FontSize',FS(10),...
        'ForegroundColor',[1 1 1],...
        'BackgroundColor',[179 179 179]/255);
    % X
    uicontrol('Parent',hReg,'Style','Text','String','x =',...
        'Position',[015 010 024 018].*WS,...
        'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
        'HorizontalAlignment','Center');
    hX = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(1)),...
        'ToolTipString','enter x-coordinate',...
        'Position',[039 010 056 020].*WS,...
        'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
        'HorizontalAlignment','Right',...
        'Tag','hX',...
        'Callback','MA_inspect_GoF(''EdWidCB'')');
    % Y
    uicontrol('Parent',hReg,'Style','Text','String','y =',...
        'Position',[100 010 024 018].*WS,...
        'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
        'HorizontalAlignment','Center')
    hY = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(2)),...
        'ToolTipString','enter y-coordinate',...
        'Position',[124 010 056 020].*WS,...
        'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
        'HorizontalAlignment','Right',...
        'Tag','hY',...
        'Callback','MA_inspect_GoF(''EdWidCB'')');
    % Z
    uicontrol('Parent',hReg,'Style','Text','String','z =',...
        'Position',[185 010 024 018].*WS,...
        'FontName',PF.times,'FontSize',FS(10),'FontAngle','Italic',...
        'HorizontalAlignment','Center')
    hZ = uicontrol('Parent',hReg,'Style','Edit','String',sprintf('%.2f',xyz(3)),...
        'ToolTipString','enter z-coordinate',...
        'Position',[209 010 056 020].*WS,...
        'FontSize',FS(10),'BackGroundColor',[.8,.8,1],...
        'HorizontalAlignment','Right',...
        'Tag','hZ',...
        'Callback','MA_inspect_GoF(''EdWidCB'')');

    % Create statistic field
    %---------------------------------------------------------------------%
    hPan = uipanel('Parent',hReg,'Title','statistic','Units','Pixels',...
        'Position',[275 005 110 040].*WS,...
        'BorderType','Beveledout', ...
        'ShadowColor',[0.5 0.5 0.5],...
        'FontAngle','Italic',...
        'FontSize',FS(10),...
        'ForegroundColor',[1 1 1],...
        'BackgroundColor',[179 179 179]/255);
    hSPM = uicontrol('Parent',hPan,'Style','Text','String','',...
        'Position',[005 001 100 020].*WS,...
        'FontSize',FS(10),...
        'HorizontalAlignment','Center');

    % Store object data
    %---------------------------------------------------------------------%
    set(hFxyz,'Tag','hFxyz','UserData',struct(...
        'hReg', [],...
        'M',    M,...
        'DIM',  DIM,...
        'XYZ',  varargin{4}.XYZmm,...
        'Z',    varargin{4}.Z,...
        'hX',   hX,...
        'hY',   hY,...
        'hZ',   hZ,...
        'hSPM', hSPM,...
        'xyz',  xyz ));
    set([hX,hY,hZ],'UserData',hFxyz);
    varargout = {hFxyz};        
 
 
%=========================================================================%
% Action: 'EdWidCB'                         callback for editable widgets %
%=========================================================================%
case 'edwidcb'
% MA_inspect_GoF('EdWidCB')
    
    % Get graphic objects
    %---------------------------------------------------------------------%
    hC    = gcbo;
    d     = find(strcmp(get(hC,'Tag'),{'hX','hY','hZ'}));
    hFxyz = get(hC,'UserData');
    UD    = get(hFxyz,'UserData');
    xyz   = UD.xyz;
    nxyz  = xyz;

    % Evaluate the ordinate
    %---------------------------------------------------------------------%
    o = evalin('base',['[',get(hC,'String'),']'],'sprintf(''error'')');
    if ischar(o) || length(o) > 1
        warning('%s: Error evaluating ordinate:\n\t%s',mfilename,lasterr);
    else
        nxyz(d) = o;
        nxyz    = spm_XYZreg('RoundCoords',nxyz,UD.M,UD.DIM);
    end

    % Check for coordinate change
    %---------------------------------------------------------------------%
    if abs(xyz(d)-nxyz(d)) > 0
        UD.xyz = nxyz;
        set(hFxyz,'UserData',UD);
        if ~isempty(UD.hReg)
            spm_XYZreg('SetCoords',nxyz,UD.hReg,hFxyz);
        end
        set(hC,'String',sprintf('%.3f',nxyz(d)));
        MA_inspect_GoF('UpdateSPMval',UD);
    end
    
    
%=========================================================================%
% Action: 'UpdateSPMval'           update SPM value in interactive window %
%=========================================================================%
case 'updatespmval'
% MA_inspect_GoF('UpdateSPMval',hFxyz)
% MA_inspect_GoF('UpdateSPMval',UD)
    
    % Update SPM value
    %---------------------------------------------------------------------%
    if nargin < 2, error('insufficient arguments'); end
    if isstruct(varargin{2}), UD = varargin{2}; else UD = get(varargin{2},'UserData'); end
    i = spm_XYZreg('FindXYZ',UD.xyz,UD.XYZ);
    if isempty(i), str = ''; else str = sprintf('%6.2f',UD.Z(i)); end
    set(UD.hSPM,'String',str);

        
%=========================================================================%
% Action: 'UpdateDataPlot'            update data plot in graphics window %
%=========================================================================%
case 'updatedataplot'
% MA_inspect_GoF('UpdateDataPlot',xyz)

    % Find co-ordinates
    %---------------------------------------------------------------------%
     xyz      = varargin{2};
    [xyz,i,c] = spm_XYZreg('NearestXYZ',xyz,GLM.XYZ);
    % dist  = sqrt(sum((GLM.XYZ - repmat(xyz,[1 GLM.w])).^2));
    % [c,i] = min(dist);        % minimal distance to current co-ordinates
    ind   = find(GLM.m==i);     % in-mask voxel index current co-ordinates
    KWyc  = GLM.KWYc(:,ind);    % measured signal at current co-ordinates
    KWzc  = GLM.KWZc(:,ind);    % predicted signal at current co-ordinates
    
    % Update data plot
    %---------------------------------------------------------------------%
    d  = 10;
    F  = spm_figure('FindWin','Graphics');
    FS = spm('FontSizes');
    figure(F)                               % figure
    subplot(2,1,1)                          % subplot
    cla reset
    switch SPM.PlotType
        
        % Refresh line graph
        %-----------------------------------------------------------------%
        case 'linegraph'
            hold on                         % data
            if ~isempty(ind)
                plot([1:GLM.n],KWyc, '-b', 'LineWidth', 1);
                plot([1:GLM.n],KWzc, '-r', 'LineWidth', 1);
                plot([1:GLM.n],mean(KWyc)*ones(1,GLM.n),  '--k', 'LineWidth', 1);
                plot([1:GLM.n],mean(KWyc)*zeros(1,GLM.n), '--g', 'LineWidth', 1);
            end;
            xlim([1 GLM.n]);                % limits
            if ~isempty(ind)
                max_amp = max(abs([max(KWyc) min(KWyc)]));
                y_limit = max_amp + (1/d)*range(KWyc);
                ylim([-y_limit +y_limit]);
                if mean(KWyc) < 0
                    text(((1/(2*d))*GLM.n),+(max_amp+(1/(2*d))*range(KWyc)),sprintf('sigma^2 = %1.4f, R^2 = %1.4f, R^2_adj = %1.4f,\nSNR_mf = %1.4f, SNR_mb = %1.4f', GLM.sig2(ind), GLM.R2(ind), GLM.adj_R2(ind), GLM.mf_SNR(ind), GLM.mb_SNR(ind)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Cap');
                    legend('measured signal', 'predicted signal', 'average signal', 'implicit baseline', 'Location', 'NorthEast');
                end;
                if mean(KWyc) > 0
                    text(((1/(2*d))*GLM.n),-(max_amp+(1/(2*d))*range(KWyc)),sprintf('SNR_mf = %1.4f, SNR_mb = %1.4f,\nsigma^2 = %1.4f, R^2 = %1.4f, R^2_adj = %1.4f', GLM.mf_SNR(ind), GLM.mb_SNR(ind), GLM.sig2(ind), GLM.R2(ind), GLM.adj_R2(ind)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');
                    legend('measured signal', 'predicted signal', 'average signal', 'implicit baseline', 'Location', 'SouthEast');
                end;
            else
                text((1/2)*GLM.n,0,'This is not an in-mask voxel!', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'FontSize', FS(24));
                ylim([-1 +1]);
            end;                            % labels
            xlabel('image', 'FontSize', FS(12));
            ylabel('value', 'FontSize', FS(12));
            title('Goodness of Fit: line graph', 'FontSize', FS(20));
            
        % Refresh scatter plot
        %-----------------------------------------------------------------%
        case 'scatterplot'
            hold on                         % data
            if ~isempty(ind)
                plot(KWzc,KWyc, 'xb', 'LineWidth', 1, 'MarkerSize', 10);
                plot([min(min([KWyc KWzc])) max(max([KWyc KWzc]))],[min(min([KWyc KWzc])) max(max([KWyc KWzc]))], '-k', 'LineWidth', 1);
            end;    
            if ~isempty(ind)                % limits
                x_range = max(KWzc) - min(KWzc);
                y_range = max(KWzc) - min(KWzc);
                xlim([min(KWzc)-(1/d)*x_range max(KWzc)+(1/d)*x_range]);
                ylim([min(KWyc)-(1/d)*y_range max(KWyc)+(1/d)*y_range]);
                axis square;
                text(min(KWzc)-(1/d)*x_range,max(KWyc)+(1/d)*y_range,sprintf('  sigma^2 = %1.4f, R^2 = %1.4f, R^2_adj = %1.4f,\n  SNR_mf = %1.4f, SNR_mb = %1.4f', GLM.sig2(ind), GLM.R2(ind), GLM.adj_R2(ind), GLM.mf_SNR(ind), GLM.mb_SNR(ind)), 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Cap');
                legend('signal plot', 'identity line', 'Location', 'SouthEast');
            else
                text((1/2)*GLM.n,0,'This is not an in-mask voxel!', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'FontSize', FS(24));
                xlim([1 GLM.n]);
                ylim([-1 +1]);
            end;                            % labels
            xlabel('predicted signal', 'FontSize', FS(12));
            ylabel('measured signal', 'FontSize', FS(12));
            title('Goodness of Fit: scatter plot', 'FontSize', FS(20));            

    end;


%=========================================================================%
% Action: 'GetCoords'                    get co-ordinates from XYZ widget %
%=========================================================================%
case 'getcoords'
% [xyz] = MA_inspect_GoF('GetCoords',hFxyz)

    % Get co-ordinates
    %---------------------------------------------------------------------%
    if nargin < 2, hFxyz = 'Interactive'; else hFxyz = varargin{2}; end
    hFxyz     = MA_inspect_GoF('FindXYZframe',hFxyz);
    varargout = {getfield(get(hFxyz,'UserData'),'xyz')};
 
 
%=========================================================================%
% Action: 'SetCoords'                      set co-ordinates to XYZ widget %
%=========================================================================%
case 'setcoords'
% [xyz,d] = MA_inspect_GoF('SetCoords',xyz,hFxyz,hC)
    
    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 4, hC = 0; else hC = varargin{4}; end
    if nargin < 3, hFxyz = MA_inspect_GoF('FindXYZframe'); else hFxyz = varargin{3}; end
    if nargin < 2, error('Set co-ords to what?'); else xyz = varargin{2}; end

    % If this is an internal call, then don't do anything
    %---------------------------------------------------------------------%
    if hFxyz == hC; return; end
    UD = get(hFxyz,'UserData');

    % Check validity of coords when called without a caller handle
    %---------------------------------------------------------------------%
    if hC <= 0
        [xyz,d] = spm_XYZreg('RoundCoords',xyz,UD.M,UD.DIM);
        if d > 0 && nargout < 2
            warning('%s: Co-ords rounded to nearest voxel centre: Discrepancy %.2f',mfilename,d);
        end
    else
        d = [];
    end

    % Update co-ordinate information & widget strings
    %---------------------------------------------------------------------%
    UD.xyz = xyz; set(hFxyz,'UserData',UD)
    set(UD.hX,'String',sprintf('%.2f',xyz(1)))
    set(UD.hY,'String',sprintf('%.2f',xyz(2)))
    set(UD.hZ,'String',sprintf('%.2f',xyz(3)))
    MA_inspect_GoF('UpdateSPMval',UD)
    MA_inspect_GoF('UpdateDataPlot',xyz)

    % Tell the registry, if we've not been called by the registry
    %---------------------------------------------------------------------%
    if (~isempty(UD.hReg) && UD.hReg ~= hC)
        spm_XYZreg('SetCoords',xyz,UD.hReg,hFxyz);
    end

    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {xyz,d};


%=========================================================================%
% Action: 'FindXYZframe'                                 find hFxyz frame %
%=========================================================================%
case 'findxyzframe'
% hFxyz = MA_inspect_GoF('FindXYZframe',h)

    % Sort out hFxyz handles
    %---------------------------------------------------------------------%
    if nargin < 2, h = 'Interactive'; else h = varargin{2}; end
    if ischar(h),  h = spm_figure('FindWin',h); end
    if ~ishandle(h),  error('invalid handle'); end
    if ~strcmp(get(h,'Tag'),'hFxyz'), h = findobj(h,'Tag','hFxyz'); end
    if isempty(h),    error('XYZ frame not found'); end
    if length(h) > 1, error('Multiple XYZ frames found'); end
    varargout = {h};


%=========================================================================%
% Action: 'Clear'                                 clear visualize figures %
%=========================================================================%
case 'clear'
% Fgraph = MA_inspect_GoF('Clear',F,mode)
% mode: 0 = clear & hide; 1 [default] = usual; 2 = RNP (?)
    
    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 3, mode = 1; else mode = varargin{3}; end
    if nargin < 2, F = 'Graphics'; else F = varargin{2}; end
    F = spm_figure('FindWin',F);

    % Clear input objects from interactive window
    %---------------------------------------------------------------------%
    % spm_input('!DeleteInputObj');

    % Get object handles from graphics window
    %---------------------------------------------------------------------%
    h = get(F,'Children');                          % Get window contents
    h = findobj(h,'flat','HandleVisibility','on');  % Drop GUI components
    g = findobj(h,'flat','Tag','PermRes');          % Look for 'PermRes'

    % If 'PermRes' object was found
    %---------------------------------------------------------------------%
    if ~isempty(g)
        tmp = get(g,'UserData');    % This has handles of permanent
        HR  = tmp.h;                % results objects in it's UserData
        HRv = tmp.g;
    else
        HR  = [];                   % No trace of permanent results objects
        HRv = {};
    end
    h = setdiff(h,HR);              % Drop permanent results objects

    % Delete stuff as appropriate
    %---------------------------------------------------------------------%
    if mode == 2                    % Don't delete axes with NextPlot
        h = setdiff(h,findobj(h,'flat','Type','axes','NextPlot','add'));
    end
    delete(h)

    % Hide the permanent stuff
    %---------------------------------------------------------------------%
    if mode == 0
        set(HR,'Visible','off')
    else
        set(HR,{'Visible'},HRv)
    end
    
    
%=========================================================================%
% Action: 'Close'                                     close visualization %
%=========================================================================%
case 'close'
% MA_inspect_GoF('Close')

    % Close visualization
    %---------------------------------------------------------------------%
    spm_clf('Interactive');
    spm_clf('Graphics');
    close(spm_figure('FindWin','Satellite'));
    % evalin('base','clear');
    cd(orig_dir);

 
%=========================================================================%
% Action: other                                     unknown action string %
%=========================================================================%
otherwise

    % Unknown action string
    %---------------------------------------------------------------------%
    error('Unknown action string!');

end