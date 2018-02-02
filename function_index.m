function function_index
% _
% List all Functions from the MACS Toolbox
% FORMAT function_index
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 15/08/2017, 17:15 (V1.1/V17)
%  Last edit: 02/02/2018, 02:55 (V1.2/V18)


% Get toolbox file index
%-------------------------------------------------------------------------%
SPM_dir  = fullfile(spm('dir'),'toolbox');
dirs     = dir(strcat(SPM_dir,'/MACS*'));
MACS_dir = fullfile(spm('dir'),'toolbox',dirs(1).name);
files    = dir(strcat(MACS_dir,'/'));

% Prepare function index
%-------------------------------------------------------------------------%
n = 0;
f = cell(1,6);

% Browse toolbox functions
%-------------------------------------------------------------------------%
for i = 1:numel(files)
    if ~isempty(strfind(files(i).name,'.m')) || ~isempty(strfind(files(i).name,'.py'))
        
        % Get function type
        %-----------------------------------------------------------------%
        n = n + 1;
        name = files(i).name(1:end-2);                  % function name
        if strncmp(name,'batch_',6), type = 'batch function';     end;
        if strncmp(name,'MA_',3),    type = 'model assessment';   end;
        if strncmp(name,'MC_',3),    type = 'model comparison';   end;
        if strncmp(name,'MS_',3),    type = 'model selection';    end;
        if strncmp(name,'ME_',3),    type = 'model estimation';   end;
        if strncmp(name,'MD_',3),    type = 'many distributions'; end;
        if strncmp(name,'MF_',3),    type = 'more functions';     end;
        if strncmp(name,'fun',3),    type = 'meta function';      end;
        
        % Read function code
        %-----------------------------------------------------------------%
        lines   = textread(strcat(MACS_dir,'/',files(i).name),'%s','delimiter','\n');
        purpose = lines{3}(3:end);                      % function purpose
        nocl    = numel(lines);                         % number of code lines
        for j = 1:numel(lines)
            if strncmp(lines{j},'% First edit:',13)     % first edit
                first_edit = lines{j}(15:end);
            end;
            if strncmp(lines{j},'%  Last edit:',13)     % last edit
                last_edit = lines{j}(15:end);
            end;
        end;
        
        % Assemble function info
        %-----------------------------------------------------------------%
        if ~isempty(strfind(files(i).name,'.m'))
            f{n,1} = name;
        else
            f{n,1} = files(i).name;
        end;
        f{n,2} = type;
        f{n,3} = purpose;
        f{n,4} = first_edit;
        f{n,5} = last_edit;
        f{n,6} = nocl;
        
    end;
end;

% Save function index
%-------------------------------------------------------------------------%
f{n+1,6} = sum(cell2mat(f(:,6)));
f{n+1,3} = 'Total Number of Code Lines';
MACS_xls = strcat(strcat(MACS_dir,'/','function_index.xls'));
xlswrite(MACS_xls, f, 'Functions', strcat('A2:F',num2str(n+2)));