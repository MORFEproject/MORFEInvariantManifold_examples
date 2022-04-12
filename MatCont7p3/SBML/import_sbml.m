function import_sbml( input_args )
%IMPORT_SBML handels gui part of importing SBML xml files
% 
global path_sys oldgds gds MC

% if sbml support not available, end here
if ~init_sbmltoolbox()
    msg = sprintf(['Unable to initialize SBML toolbox\n'...
           'See command window for more information']);
    errordlg(msg,'Importing SBML file', 'modal');
    return;
end

gds = close_windows(MC, gds);

% set filter_spec/default_sbml_path to a sane default
% get /SBML dir
sbml_dir = mfilename('fullpath');
sbml_dir = sbml_dir(1:end-length(mfilename));

sbml_settings_file = fullfile(sbml_dir,'settings.mat');
if exist(sbml_settings_file,'file')
    load(sbml_settings_file,'default_sbml_path')
end
if ~exist('default_sbml_path','var') || ~exist(default_sbml_path,'dir')
    default_sbml_path = '.';
end
filter_spec = fullfile(default_sbml_path,'*.xml');

% ask file to user
[file_name, path_name] = uigetfile(filter_spec,'Select SBML file');

% if cancel selected, end here
if isnumeric(file_name)
    open_windows(MC, gds);
    return
end

% update sane filter_spec/default_sbml_path
default_sbml_path = path_name;
if exist(sbml_settings_file,'file')
    save(sbml_settings_file,'default_sbml_path','-append');
else
    save(sbml_settings_file,'default_sbml_path');
end

% get sbml_model from file
sbml_file = fullfile(path_name,file_name);
model = TranslateSBML(sbml_file); % catch error?


% show dialog that we are importing the file
% if ~isempty(whos('global','NDEBUG')) % to use, declare 'global NDEBUG'
%     dlg = dialog('Name', 'Importing SBML file', 'WindowStyle', 'normal');
% else
%     dlg = dialog('Name', 'Importing SBML file', 'WindowStyle', 'modal');
% end
% dlg_pos = get(dlg,'Position');
% dlg_pos(3:4) = [300 75];
% set(dlg, 'Position', dlg_pos);
% interior_pos = get(dlg,'Position');
% interior_pos(1:2) = [20 20];
% interior_pos(3:4) = interior_pos(3:4) - 40;
% dlg_message = sprintf('Please wait, lengthy operation...');
% uicontrol(dlg, 'Position', interior_pos, 'Style', 'edit', ...
%     'Max', 5, 'String', dlg_message);
% drawnow

h = waitbar(0,'Importing SBML file');

% get the gds structure from sbml, might be a lengthy operation
if ~isempty(whos('global','NDEBUG'))
    gds_struct = gds_from_sbml(model);
else
    try
        gds_struct = gds_from_sbml(model);
    catch exception
        delete(h);
        errordlg(exception.message, ...
            'Unable to convert SBML file to matcont', 'modal');
        open_windows(MC, gds);
        return;
    end
end

% hard work is done, remove modal dialog from screen
%delete(dlg);
delete(h);

% store the oldgds, enable cancel to work in system dialog
if isfield(gds,'system') && ~isempty(gds.system)
    oldgds = gds;
end

% set global gds structure
gds = gds_struct;

system_mat_file = fullfile(path_sys,gds.system);
save(system_mat_file,'gds');
systems;

end

function gds = close_windows(MC, gds)
% CLOSE_WINDOWS closes all open windows except the main matcont window
%               and updates the gds structure
% TODO exactly the same code is used in systems.m and editsystem.m
%      but a difference between system and editsystem, what is right way
    %starter-window open?
    if ~isempty(MC.starter), gds.open.figuur = 1;else gds.open.figuur = 0;end
    close(MC.starter);
    % continuer-window open?
    if ~isempty(MC.continuer), gds.open.continuer = 1;else gds.open.continuer = 0;end
    close(MC.continuer);
    % numeric window open?
    if ~isempty(MC.numeric_fig), gds.open.numeric_fig = 1;else gds.open.numeric_fig = 0;end
    close(MC.numeric_fig);
    %2D-plot open
    if ~isempty(MC.D2), gds.open.D2 = size(MC.D2);else gds.open.D2 = 0;end
    close(MC.D2);
    %3D-plot open
    if ~isempty(MC.D3), gds.open.D3 = size(MC.D3);else gds.open.D3 = 0;end
    close(MC.D3);
    %3D-plot open
    if ~isempty(MC.PRC), gds.open.PRC = size(MC.PRC);else gds.open.PRC = 0;end
    try
        close(MC.PRC);
        gds.open.PRC = 0;
    catch
        MC.PRC = [];
    end
    %3D-plot open
    if ~isempty(MC.dPRC), gds.open.dPRC = size(MC.dPRC);else gds.open.dPRC = 0;end
    try
        close(MC.dPRC);
        gds.open.dPRC = 0;
    catch
        MC.dPRC = [];
    end
    if ~isempty(MC.integrator),gds.open.integrator = 1;else gds.open.integrator = 0;end;
    close(MC.integrator);
end

function open_windows(MC, gds)
% OPEN_WINDOWS opens all windows that gds says should be open
% TODO exactly the same code is used in systems.m and editsystem.m
    if gds.open.figuur==1;starter;end
    if gds.open.continuer==1;continuer;end
    if gds.open.numeric_fig==1;numeric;end
    if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
    if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
    if gds.open.integrator==1;integrator;end
end
