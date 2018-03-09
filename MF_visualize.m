function varargout = MF_visualize(varargin)
% _
% Visualize High-Dimensional Data using SPM Overlay Mode
% FORMAT varargout = MF_visualize(varargin)
% 
%     varargin  - cell array of input variables
%     varargout - cell array of output variables
% 
% FORMAT xVis = MF_visualize('Specify') opens the interactive mode to
% specify a visualization scheme including the data to be visualized, the
% desired overlay image, an overlay threshold and some plotting parameters.
% Together, this is returned as the persistent structure variable xVis.
% 
%     xVis - structure with visualization parameters; required fields:
%     o data     - N x 1 cell array with images pathes to visualize
%     o overlay  - 1 x 1 cell array with image path to overlay image
%     o thresh   - string indicating threshold (e.g. '>0.5' or '==64')
%     o PlotType - string indicating plot type ('bar', 'plot' or 'matrix')
%     o LineSpec - string with specification according to 'doc LineSpec'
%     o XTicks   - 1 x N cell array of strings with labels for the x-axis
%     o YLimits  - 1 x 2 vector of doubles with limits for the y-axis
%     o Title    - string, title of the plot
% 
% FORMAT xDat = MF_visualize('LoadData',xVis) uses a visualization scheme
% xVis and loads input images, overlay image and some technical parameters.
% Together, this is returned as the persistent structure variable xDat.
% 
%     xDat - structure with information on data; required fields:
%     o D   - N x V data matrix with visualized data
%     o N   - integer, number of images
%     o V   - integer, number of voxels
%     o H   - structure with overlay header
%     o O   - 1 x V vector with overlay image data
%     o XYZ - 3 x V vector with overlay co-ordinates
% 
% FORMAT [hReg,xSPM,SPM] = MF_visualize('Setup',xVis) uses a visualization
% scheme xVis, loads the corresponding data xDat and prepares figures,
% graphics and handles for data visualization. The function should be
% called in this mode when xVis is specified externally and not using the
% SPM interface.
% 
% This function was written by modifying the SPM Results User Interface.
% Consequently, much of this code is stolen from spm_results_ui.m!
% Additional inspiration came from spm_run_bms_vis.m.
% 
% References:
% [1] Soch J, Allefeld C (2018): "MACS - a new SPM toolbox for model
%     assessment, comparison and selection". Journal of Neuroscience
%     Methods, in review. URL: https://www.biorxiv.org/content/early/2017/11/09/194365.
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
persistent xVis xDat orig_dir
% xVis     - visualization scheme  (see comments above)
% xDat     - data to be visualized (see comments above)
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
% xVis = MF_visualize('Specify')

    % Initialise
    %---------------------------------------------------------------------%
    Finter = spm('FigName','MF_visualize: Specify');
    
    % Select images and overlay
    %---------------------------------------------------------------------%
    xVis.data    = cellstr(spm_select(Inf, 'image', 'Please select your input images!', [], pwd, '.*', '1'));
    xVis.overlay = cellstr(spm_select(1, 'image', 'Please select the overlay image!', [], pwd, '.*', '1'));
    xVis.thresh  = spm_input('Overlay threshold:',1,'s','>0.5');
    
    % Enter plotting parameters
    %---------------------------------------------------------------------%
    PlotTypes     = {'bar','plot','matrix'};
    xVis.PlotType = PlotTypes{spm_input('Plot Type:','+1','b',{'bar','plot','matrix'},[1 2 3])};
    xVis.LineSpec = spm_input('Line Spec:','+1','s','b');
    xVis.XTicks   = spm_input('X-Axis Ticks:','+1','e','{}');
    xVis.YLimits  = spm_input('Y-Axis Limits:','+1','r','[]');
    xVis.Title    = spm_input('Plot Title:','+1','s','Title');
    clear PlotTypes

    % Set up visualization
    %---------------------------------------------------------------------%
    MF_visualize('Setup',xVis);

    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {xVis};
    

%=========================================================================%
% Action: 'Setup'                                    set up visualization %
%=========================================================================%
case 'setup'
% [hReg,xSPM,SPM] = MF_visualize('Setup',xVis)

    % Return
    %---------------------------------------------------------------------%
    if nargin < 2
        xVis = MF_visualize('Specify');
        return
    end;
    
    % Initialise
    %---------------------------------------------------------------------%
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MF_visualize: Setup');
    spm_clf('Satellite');
    orig_dir = pwd;
 
    % Load input images
    %---------------------------------------------------------------------%
    xVis = varargin{2};
    xDat = MF_visualize('LoadData',xVis);

    % Find activated voxels
    %---------------------------------------------------------------------%
    eval(strcat('stv_ind = find(xDat.O',xVis.thresh,');')); 
    XYZ_ord      = zeros(3,xDat.V);         % ordinal voxel co-ordinates
    XYZ_ord(1,:) = kron(ones(1,xDat.H.dim(3)*xDat.H.dim(2)), [1:xDat.H.dim(1)]);
    XYZ_ord(2,:) = kron(ones(1,xDat.H.dim(3)), kron([1:xDat.H.dim(2)], ones(1,xDat.H.dim(1))));
    XYZ_ord(3,:) = kron([1:xDat.H.dim(3)], ones(1,xDat.H.dim(2)*xDat.H.dim(1)));
    XYZ_ind      = XYZ_ord(:,stv_ind);      % supra-threshold voxel indexes
    
    % Create xSPM structure
    %---------------------------------------------------------------------%
    xSPM.swd       = fileparts(xVis.overlay{1});
    xSPM.title     = xVis.Title;
    xSPM.Z         = xDat.O(stv_ind);
    xSPM.M         = xDat.H.mat;
    xSPM.iM        = inv(xSPM.M);
    xSPM.VOX       = sqrt(sum(xSPM.M(1:3,1:3).^2));
    xSPM.DIM       = xDat.H.dim(1:3)'; 
    xSPM.XYZ       = XYZ_ind;
    xSPM.XYZmm     = xSPM.M * [xSPM.XYZ; ones(1,size(xSPM.XYZ,2))];
    xSPM.XYZmm     = xSPM.XYZmm(1:3,:);
    xSPM.units     = {'mm' 'mm' 'mm'};
    xSPM.thresDesc = 'none';
    
    % Create SPM structure
    %---------------------------------------------------------------------%
    SPM.swd = xSPM.swd;
    cd(SPM.swd);
 
    % Get space information
    %---------------------------------------------------------------------%
    M   = xSPM.M;
    DIM = xSPM.DIM;
    
    % Set up visualize GUI
    %---------------------------------------------------------------------%
    spm_clf(Finter);
    spm('FigName','MF_visualize: Setup',Finter);
    hReg = MF_visualize('SetupGUI',M,DIM,xSPM,Finter);

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
    varargout = {hReg,xSPM,SPM};
    spm('Pointer','Arrow');


%=========================================================================%
% Action: 'LoadData'                                    load input images %
%=========================================================================%
case 'loaddata'
% xDat = MF_visualize('LoadData',xVis)
    
    % Get data dimensions
    %---------------------------------------------------------------------%
    xVis   = varargin{2};                   % xVis structure
    xDat.H = spm_vol(xVis.overlay{1});      % overlay header
    xDat.V = prod(xDat.H.dim);              % number of voxels
    xDat.N = size(xVis.data,1);             % number of volumes
    
    % Load input images
    %---------------------------------------------------------------------%
    spm_progress_bar('Init', 100, 'Load data...' , '');
    xDat.D = zeros(xDat.N,xDat.V);          % data matrix
    for i = 1:xDat.N
        D_hdr = spm_vol(xVis.data{i});      % data header
        D_img = spm_read_vols(D_hdr);       % data image
        xDat.D(i,:) = reshape(D_img,1,[]);  % data vector
        spm_progress_bar('Set',(i/xDat.N)*100);
    end;
    clear D_hdr D_img
    spm_progress_bar('Clear');
    
    % Load overlay image
    %---------------------------------------------------------------------%
    [O_img XYZ] = spm_read_vols(xDat.H);
    xDat.O   = reshape(O_img,[1 xDat.V]);   % overlay image
    xDat.XYZ = XYZ;                         % XYZ co-ordinates
    clear O_img XYZ
    
    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {xDat};


%=========================================================================%
% Action: 'SetupGUI'                                 set up visualize GUI %
%=========================================================================%
case 'setupgui'
% hReg = MF_visualize('SetupGUI',M,DIM,xSPM,Finter)

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
    hFxyz      = MF_visualize('DrawXYZgui',M,DIM,varargin{4},xyz,hReg);
    spm_XYZreg('XReg',hReg,hFxyz,'MF_visualize');

    % Setup buttons for visualization functions
    %---------------------------------------------------------------------%
    MF_visualize('DrawButts',hReg,DIM,Finter,WS,FS);
    set(findobj(hPan),'Units','Normalized','FontUnits','Normalized');

    % Return output variables
    %---------------------------------------------------------------------%
    varargout = {hReg};
 

%=========================================================================%
% Action: 'DrawButts'                  draw buttons in interactive window %
%=========================================================================%
case 'drawbutts'
% MF_visualize('DrawButts',hReg,DIM,Finter,WS,FS)

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
        'Callback',['MF_visualize(''Clear''); ',...
          'spm_input(''!DeleteInputObj''),',...
          'spm_clf(''Satellite'')'],...
        'Interruptible','on','Enable','on',...
        'DeleteFcn','spm_clf(''Graphics'')',...
        'Position',[280 050 048 020].*WS);
    uicontrol('Parent',hReg,'Style','PushButton','String','exit',...
        'ToolTipString','Exit the results section',...
        'FontSize',FS(9),'ForegroundColor','r',...
        'Callback','MF_visualize(''close'')',...
        'Interruptible','on','Enable','on',...
        'Position',[332 050 048 020].*WS);

    
%=========================================================================%
% Action: 'DrawXYZgui'                 draw XYZ GUI in interactive window %
%=========================================================================%
case 'drawxyzgui'
% hFxyz = MF_visualize('DrawXYZgui',M,DIM,xSPM,xyz,hReg)

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
        'Callback','MF_visualize(''EdWidCB'')');
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
        'Callback','MF_visualize(''EdWidCB'')');
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
        'Callback','MF_visualize(''EdWidCB'')');

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
% MF_visualize('EdWidCB')
    
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
        MF_visualize('UpdateSPMval',UD);
    end
    
    
%=========================================================================%
% Action: 'UpdateSPMval'           update SPM value in interactive window %
%=========================================================================%
case 'updatespmval'
% MF_visualize('UpdateSPMval',hFxyz)
% MF_visualize('UpdateSPMval',UD)
    
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
% MF_visualize('UpdateDataPlot',xyz)

    % Find co-ordinates
    %---------------------------------------------------------------------%
     xyz      = varargin{2};
    [xyz,i,c] = spm_XYZreg('NearestXYZ',xyz,xDat.XYZ);
    % dist  = sqrt(sum((xDat.XYZ - repmat(xyz,[1 xDat.V])).^2));
    % [c,i] = min(dist);        % minimal distance to current co-ordinates
    
    % Update data plot
    %---------------------------------------------------------------------%
    F  = spm_figure('FindWin','Graphics');
    FS = spm('FontSizes');
    figure(F)                               % figure
    subplot(2,1,1)                          % subplot
    switch xVis.PlotType
        
        % Refresh bar plot
        %-----------------------------------------------------------------%
        case 'bar'
            bar(xDat.D(:,i),xVis.LineSpec);
            xlim([(1-0.5-2*(1/xDat.N)) (xDat.N+0.5+2*(1/xDat.N))]);
            if ~isempty(xVis.XTicks)        % XTicks
                set(gca,'XTick',[1:xDat.N],'XTickLabel',xVis.XTicks);
            else
                set(gca,'XTick',[1:xDat.N]);
            end;
        
        % Refresh line plot
        %-----------------------------------------------------------------%
        case 'plot'
            plot([1:xDat.N],xDat.D(:,i),xVis.LineSpec,'LineWidth',1);
            xlim([1 xDat.N]);
        
        % Refresh matrix plot
        %-----------------------------------------------------------------%
        case 'matrix'
            imagesc(reshape(xDat.D(:,i),[sqrt(xDat.N) sqrt(xDat.N)]));
            axis([(1-0.5) (sqrt(xDat.N)+0.5) (1-0.5) (sqrt(xDat.N)+0.5)]);
            axis ij
            axis square
            colorbar('Location','EastOutside');
            if ~isempty(xVis.XTicks)        % XTicks
                set(gca,'XTick',[1:sqrt(xDat.N)],'XTickLabel',xVis.XTicks);
                set(gca,'YTick',[1:sqrt(xDat.N)],'YTickLabel',xVis.XTicks);
            else
                set(gca,'XTick',[1:sqrt(xDat.N)]);
                set(gca,'YTick',[1:sqrt(xDat.N)]);
            end;
        
    end;
    
    % Label data plot
    %---------------------------------------------------------------------%
    if ~isempty(xVis.YLimits)               % YLimits
        ylim(xVis.YLimits);
    end;
    if ~strcmp(xVis.PlotType,'matrix')
        xlabel('image', 'FontSize', FS(12));% labels
        ylabel('value', 'FontSize', FS(12));
    end;
    title(xVis.Title, 'FontSize', FS(20));  % title


%=========================================================================%
% Action: 'GetCoords'                    get co-ordinates from XYZ widget %
%=========================================================================%
case 'getcoords'
% [xyz] = MF_visualize('GetCoords',hFxyz)

    % Get co-ordinates
    %---------------------------------------------------------------------%
    if nargin < 2, hFxyz = 'Interactive'; else hFxyz = varargin{2}; end
    hFxyz     = MF_visualize('FindXYZframe',hFxyz);
    varargout = {getfield(get(hFxyz,'UserData'),'xyz')};
 
 
%=========================================================================%
% Action: 'SetCoords'                      set co-ordinates to XYZ widget %
%=========================================================================%
case 'setcoords'
% [xyz,d] = MF_visualize('SetCoords',xyz,hFxyz,hC)
    
    % Retrieve input arguments
    %---------------------------------------------------------------------%
    if nargin < 4, hC = 0; else hC = varargin{4}; end
    if nargin < 3, hFxyz = MF_visualize('FindXYZframe'); else hFxyz = varargin{3}; end
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
    MF_visualize('UpdateSPMval',UD)
    MF_visualize('UpdateDataPlot',xyz)

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
% hFxyz = MF_visualize('FindXYZframe',h)

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
% Fgraph = MF_visualize('Clear',F,mode)
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
% MF_visualize('Close')

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