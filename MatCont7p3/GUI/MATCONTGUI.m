function GUI(x)
% This function contructs the MatCont GUI. It creates 'Session' and launches the Main Window. This function also creates the menu of the main window.
% x is a number and is the number of arguments passed along the main driver function.
%
% x == 0: no arguments -> restore state
% x == 1: one argument -> do not restore state
% x == 2: two arguments -> clear all and do not restore state (used to clear cached code from memory)
%.


global session;
if x == 2
    clear all;  %clear all, including cached code
    delete(allchild(0)); %this kills all MATLAB windows
    x = 2;
end


%Construction of the Main model Object. This object is never swapped, only reorganizes.
session = Session();

if x == 0
    % no arguments given, restoring state...

    % './Systems/session.mat' contains the previous state of session.
    % Session will now load in the state and return two values, one to restore the open windows and one to restore plotwindows.
    % this step has to occur after the GUI Mainwindow has been constructed.
    [restorewindows, plotconfs] = session.loadFromFile();
else
   % load in empty dummies if no restoration needs to happen
   restorewindows = @(s)[];
   plotconfs = {};
end



%Debug only, enable to monitor Session-events
%{
session.addlistener('settingChanged', @(o, e) fprintf('*** SESSION: settingChanged\n'));
session.addlistener('settingsChanged', @(o, e) fprintf('*** SESSION: settingsChanged\n'));
session.addlistener('computationChanged', @(o, e) fprintf('*** SESSION: computationChanged\n'));
session.addlistener('solutionChanged', @(o, e) fprintf('*** SESSION: solutionChanged\n'));
session.addlistener('initpointChanged', @(o, e) fprintf('*** SESSION: initpointChanged\n'));
session.addlistener('lockChanged', @(o, e) fprintf('*** SESSION: lockChanged\n'));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% Starting to build the main window %% %%
windowmanager = session.getWindowManager();
mwindowhandle = windowmanager.createMainWindow();


handlelist_computing = []; %enter handles that lock when computing.
handelist_systemneeded = []; %enter handles that disable when no system is selected

% % % Select MENU % % %
handle = uimenu(mwindowhandle , 'Label' , 'Select');
handlelist_computing(end+1) = handle;

syshandle = uimenu(handle , 'Label' , 'System');
handlelist_computing(end+1) = syshandle;

GUIMenuItem(syshandle, @(o,e) GUIMainPanel.call_new_system(session), @() 'New',...
    [] , @() session.isUnlocked() , '');

GUIMenuItem(syshandle ,@(o,e) CLBrowsers.systemBrowser(session) , @() 'Load/Edit/Delete Systems',...
    [] , @() session.isUnlocked() , '', 'Separator' , 'on');

GUIMenuItem(syshandle , @(o,e) GUIMainPanel.call_userfunctions(session) , @() 'Manage Userfunctions' , session,...
                @() session.hasSystem() && session.isUnlocked() , 'settingsChanged');
          
            
GUIMenuItem(handle ,@(o,e) CLBrowsers.diagramBrowser(session), @() 'Diagram',...
    session , @() session.hasSystem() && session.isUnlocked() , 'settingsChanged', 'Separator' , 'on' );
            
GUIMenuItem(handle ,@(o,e) CLBrowsers.solutionBrowser(session), @() 'Curve',...
    session , @() session.hasSystem() && session.isUnlocked() , 'settingsChanged', 'Separator' , 'on' );

GUIMenuItem(handle, @(o,e) CLBrowsers.switchBrowser(session), @() 'Initial Point',...
    session , @() session.hasSolution() && session.isUnlocked() , 'solutionChanged', 'Separator' , 'on');


   
GUIMenuItem(handle, @()  GUIDiagramOrganizerModel.manageCurrentDiagram(session)  , 'Organize Diagrams' , session , @() session.hasSystem() , 'settingsChanged', 'Separator' , 'on');


%add the kill-switch. GUI shuts down when the mainwindow handle receives a delete-signal.
uimenu(handle, 'Label' , 'Exit' , 'Callback' , @(o,e) delete(mwindowhandle), 'Separator' , 'on');


% % % % Type Menu % % %  
handle = uimenu(mwindowhandle , 'Label' , 'Type');
handlelist_computing(end+1) = handle;
handelist_systemneeded(end+1) = handle;

%specialized menu items: Select Point and Select Curve.
GUISelectPointsMenu(handle, session);
GUISelectCurveMenu(handle, session, session.branchmanager);



% % % Window/Output Menu % % %
handle = uimenu(mwindowhandle , 'Label' , 'Window/Output');
handlelist_computing(end+1) = handle;
handelist_systemneeded(end+1) = handle;
gm = uimenu(handle , 'Label' , 'Graphic');

createplot2 = uimenu(gm, 'Label' , '2D plot' , 'Callback' , @(o,e) session.getOutputManager().createPlot(2), 'Separator', 'on');
createplot3 = uimenu(gm, 'Label' , '3D plot' , 'Callback' , @(o,e) session.getOutputManager().createPlot(3));
GUIPreviousPlotsMenuItem(gm, 'uimenu', session, 'Separator', 'on');
GUIWindowLaunchButton(handle , windowmanager , 'uimenu', ...
       @(fhandle) session.getOutputManager().createNumeric(fhandle) , 'numeric'); 
   
GUIWindowLaunchButton(handle , session.getWindowManager() , 'uimenu', ...
       @(fhandle) GUIRenderSettings(fhandle, session, [DefaultValues.SECTIONID.CONTINUER, DefaultValues.SECTIONID.INTEGRATOR]), 'continuer','Separator', 'on' );
GUIWindowLaunchButton(handle , session.getWindowManager() , 'uimenu', ...
       @(fhandle)  GUIRenderSettings(fhandle, session, DefaultValues.SECTIONID.STARTER)  , 'starter');   
   
   
% % % Compute Menu % % %
%specialized menu: compute
menu = GUIComputeMenu(mwindowhandle, 'uimenu', session);

   


 % % % Options Menu % % %  

handle = uimenu(mwindowhandle , 'Label' , 'Options');
handlelist_computing(end+1) = handle;
handelist_systemneeded(end+1) = handle;

%Add options menu for the current 'settings' in session: 
% change here the options in the menu. Make sure these options are always available in every 'settings' (define in 'CompConf')
GUIOptionsMenu.install(handle, session, {'option_pause', 'option_archive', 'option_output'});


%handlelist_computing(end+1) = handle;

GUIWindowLaunchButton(handle , session.getWindowManager() , 'uimenu', ...
       @(fhandle) GUILineConfigPanel(fhandle , session) , 'coloropt', 'separator', 'on');   
  
%MENU version, disabled.   
% subhandle = uimenu(handle , 'Label' , 'Computational Options', 'separator', 'on');
% option has to always exist in the current 'settings' in session (define in 'CompConf')
% GUIOptionsMenu.install(subhandle, session, {'option_increment', 'option_moorepenrose', 'option_tsearchorder'});

% % % % % % % Help Menu % % % % %

handle = uimenu(mwindowhandle , 'Label' , 'Help');
uimenu(handle , 'Label' , 'MatCont Help' , 'Callback' , @(o,e) web(which('doc_matcont.html')));
uimenu(handle , 'Label' , 'Emergency Reset GUI' , 'Callback' , @(o,e) session.resetState());

%%% %%%     construct the main panel     %%%   %%%%%%%%%%%
mainpanel = GUIMainPanel(mwindowhandle, session);

%make sure the labels in the list are disabled when session is locked (computing)
GUILockEnableSync(session, handlelist_computing, 'lockChanged', @() session.isUnlocked());
%make sure the labels in the list are disabled when no system is loaded
GUILockEnableSync(session, handelist_systemneeded, 'settingsChanged', @() session.hasSystem());
GUILockEnableSync(session, [createplot2, createplot3], 'computationChanged', @() ~isempty(session.getOutputInterpreter()));

%reveal end result.
set(mwindowhandle, 'Visible', 'on');


%Add listener: open starter and continuer/integrator on 'computation' change. (if not
%already open)

session.addlistener('computationChanged', @(o,e) openStartCont(session));

%restore windows-status from previous state  (if any)
restorewindows(session);

%restore plots from previous state (if any)
session.outputmanager.loadPlotConfs(plotconfs);

%done
end



function openStartCont(session)
    %open starter/continuer/integrator if not already open
    if ~session.windowmanager.isWindowOpen('starter')
        GUIRenderSettings(session.windowmanager.demandWindow('starter')  , session, DefaultValues.SECTIONID.STARTER);
    end
    if ~session.windowmanager.isWindowOpen('continuer')
        GUIRenderSettings(session.windowmanager.demandWindow('continuer')  , session, [DefaultValues.SECTIONID.CONTINUER, DefaultValues.SECTIONID.INTEGRATOR]);
    end
end

