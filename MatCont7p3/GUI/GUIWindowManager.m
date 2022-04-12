classdef GUIWindowManager < handle
    properties
        allowedfunction
        windowhandles
        windownames
        
        enabledvalues
        
        session;
        eventlistener;
        eventlistener2;
        
        windowpositions;
    end
    events
       windowChanged 
    end
    methods
        function obj = GUIWindowManager(session)
            if ~isempty(session)
                
                obj.windowpositions = session.globalsettings.getSetting('windowpositions');
            else
                obj.windowpositions = CLSettingWindowPositions();
            end
            
            obj.allowedfunction = struct();
            obj.windowhandles = struct();
            obj.windownames = struct();
            obj.session = session;
            
            obj.registerWindow('main', 'MatCont GUI', @() 1, [-1000 -1000 640 270]);
            obj.registerWindow('numeric', 'Numeric', @() session.hasSystem(), [-1000 -1000 300   450]);
            obj.registerWindow('starter', 'Starter', @() session.hasSystem(), [-1000 -1000 320   530]);
            obj.registerWindow('continuer', 'Continuer/Integrator', @() session.hasSystem(), [-1000 -1000 380   530]);
            obj.registerWindow('monitor', 'Control', @() session.hasSystem(),  [-1000 -1000 340 116]);
            obj.registerWindow('browser', 'Data Browser', @()   session.isUnlocked(), [ -1000 -1000 680 400]);
            obj.registerWindow('plotlayout2d', 'Layout', @()   session.isUnlocked(), [-1000 -1000  400 180]);
            obj.registerWindow('plotlayout3d', 'Layout', @()  session.isUnlocked(), [-1000 -1000  400 280]);
            obj.registerWindow('diagramorganizer', 'Organize Diagram',  @() session.hasSystem() && session.isUnlocked(), [-1000 -1000  400 370]);
            obj.registerWindow('coloropt', 'Plot Properties', @() session.isUnlocked(), [ -1000 -1000 600   450]);
            
            obj.registerWindow('PRC', 'PRC', @() 1, [ -1000 -1000 600   450]);
            obj.registerWindow('dPRC', 'dPRC', @() 1, [ -1000 -1000 600   450]);
            
            if ~isempty(session)
                obj.eventlistener = session.addlistener('settingsChanged' , @(o,e) obj.sessionStateChanged()); 
                obj.eventlistener2 = session.addlistener('lockChanged' , @(o,e) obj.sessionStateChanged()); 
            end
            
        end
        
        function closeLayoutWindows(obj)
            obj.closeWindow('plotlayout2d');
            obj.closeWindow('plotlayout3d');
        end
        
        function sessionStateChanged(obj)
            names = fieldnames(obj.windowhandles);
            dirty = 0;
            for i = 1:length(names)
                newval = feval(obj.allowedfunction.(names{i}));
                if(obj.enabledvalues.(names{i}) ~= newval)
                    obj.enabledvalues.(names{i}) = newval;
                    dirty = 1;
                    if ((~newval) && obj.isWindowOpen(names{i}))
                        obj.closeWindow(names{i});
                    end
                    
                end
            end
            if dirty
                obj.notify('windowChanged');
            end
        end

        
        
        function registerWindow(obj, label, name, func, defaultposition)
            obj.allowedfunction.(label) = func;
            
            obj.windowpositions.setDefault(label, defaultposition);
            obj.windownames.(label) = name;
            obj.windowhandles.(label) = [];
            obj.enabledvalues.(label) = -1;
            
        end
        function b = isAllowed(obj, label)
            b = feval(obj.getAllowed(label));
        end
        function f = getAllowed(obj, label)
            if isfield(obj.allowedfunction, label)
                f = obj.allowedfunction.(label);
            else
                fprintf(2, 'WM: unknown label: %s\n', label);
                f = @() 1;
            end
        end
        
        function b = isWindowOpen(obj, label)
            b = ~isempty(obj.windowhandles.(label)) && isvalid(obj.windowhandles.(label));
        end
        function b = isWindowOpenable(obj, label)
            b =  ~obj.isWindowOpen(label) && feval(obj.allowedfunction.(label));
        end
        function closeWindow(obj, label)
            if obj.isWindowOpen(label)
                obj.setWindowPosition(label);
                delete(obj.windowhandles.(label));
            end
            obj.windowhandles.(label) = [];
            obj.notify('windowChanged');
        end
        
        function name = getWindowName(obj, label)
            name = obj.windownames.(label);
        end
        
        function window = createWindow(obj, label)
            window = obj.createFigure();
            set(window,'CloseRequestFcn', @(src,ev) closeSubwindow(obj, label, src) , 'DeleteFcn' , @(src,ev) deleteWindow(src,ev), 'Name' , obj.getWindowName(label));
            obj.windowhandles.(label) = window;
            obj.installPosition(label);
            obj.notify('windowChanged');
        end
        
        function f = createFigure(~)
            f = figure('Visible' , 'off' , 'NumberTitle' , 'off','MenuBar' , 'none', 'tag', 'matcont', 'KeyPressFcn', @(o,e) SessionOutputKeypress(o, e));
        end
        
        function window = demandWindow(obj, label)
            obj.closeWindow(label);
            window = obj.createWindow(label);
        end
        
        
        
        function window = createMainWindow(obj)
           window = obj.createWindow('main');
           set(window, 'CloseRequestFcn', @(src,ev) closeMainwindow(obj, src));
        end
        
        function closeSubWindows(obj)
            windowtags = fieldnames(obj.windowhandles);
            for i = 1:length(windowtags)
                if ~strcmp(windowtags{i}, 'main')
                    obj.closeWindow(windowtags{i});
                end
            end
        end
        
        function closeAllWindows(obj)
            if ~obj.isWindowOpen('main'); obj.closeSubWindows(); return; end
            mainhandle = obj.windowhandles.main;
            
            set(mainhandle, 'Visible' , 'off');
            drawnow;
            
            obj.closeSubWindows();
            obj.closeWindow('main');
        end
        
        
        function pos = getWindowPosition(obj, label)
            pos = obj.windowpositions.getPosition(label);
        end
        function setWindowPosition(obj, label)
            obj.windowpositions.setPosition(label, get(obj.windowhandles.(label), 'Position'));
            
        end
        
        function installPosition(obj, label)
            if ~obj.isWindowOpen(label); return; end
            
            whandle = obj.windowhandles.(label);
            pos = obj.getWindowPosition(label);
            
            set(whandle , 'Position' , pos);
            if (pos(1) < -100)
                movegui(whandle, 'center');
            else
                movegui(whandle, 'onscreen');
            end
        end
        
        
        function ff = restoreFunction(obj)
            ff = @(s) 1;
            if obj.isWindowOpen('starter')
                ff = @(s) {ff(s), GUIRenderSettings(s.windowmanager.demandWindow('starter')  , s, DefaultValues.SECTIONID.STARTER)};
            end
            if obj.isWindowOpen('continuer')
                ff = @(s) {ff(s), openContinuer(s, DefaultValues.SECTIONID.CONTINUER)};
                ff = @(s) {ff(s), openContinuer(s, DefaultValues.SECTIONID.INTEGRATOR)};
            end
            if obj.isWindowOpen('numeric')
                ff = @(s) {ff(s), launchNumeric(s)};
            end
        end
        
    end
    
end

function handle = openContinuer(s, sectiontoken)
if ~s.windowmanager.isWindowOpen('continuer')
    handle = GUIRenderSettings(s.windowmanager.demandWindow('continuer')  , s, sectiontoken);
else
    handle = []; 
end
end


function handle = launchNumeric(s)
    handle = s.windowmanager.demandWindow('numeric');
    s.getOutputManager().createNumeric(handle);
    set(handle, 'Visible', 'on');

end


function closeSubwindow(windowmanager, tag, src)
if isvalid(windowmanager)
    windowmanager.closeWindow(tag);
else
    delete(src);
end
end

function closeMainwindow(windowmanager, src)
if isvalid(windowmanager)
    windowmanager.closeWindow('main');
    %windowmanager.closeAllWindows();
else
   delete(src); 
end

end
function deleteWindow(src,~)
figure(src); %INPUT BUG FIX
delete(src);
end
