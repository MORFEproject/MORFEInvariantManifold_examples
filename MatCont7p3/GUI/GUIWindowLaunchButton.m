classdef GUIWindowLaunchButton < handle
    
    properties
        handle
        eventlistener
        
        windowmanager
        launchfunction
        windowtag;
        
    end
    
    methods
        function obj =  GUIWindowLaunchButton(parenthandle , wmanager,  uitype, launchfunction , windowtag , varargin)
            obj.windowmanager = wmanager;
            obj.launchfunction = launchfunction;
            obj.windowtag = windowtag;
            
            if (strcmp(uitype,'uimenu'))
               uiconstructor = @uimenu;
               vars = {'Label' , wmanager.getWindowName(windowtag)};
            else
                uiconstructor = @uicontrol;
                vars = {'String' , wmanager.getWindowName(windowtag) , 'Style' , 'pushbutton'};
            end
            
            
            obj.handle = uiconstructor( parenthandle , vars{:} , 'Callback' , @(o,e) obj.launchWindow() , varargin{:} );
            
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.eventlistener = obj.windowmanager.addlistener('windowChanged' , @(o,e) obj.windowChanged());
        
            obj.windowChanged();
        end
        
        function launchWindow(obj)
            if (obj.windowmanager.isWindowOpenable( obj.windowtag))
                f = obj.windowmanager.createWindow(obj.windowtag);
                obj.launchfunction(f);
                if ishandle(f); set(f,'Visible' , 'on'); end
            end
        end
        
        function windowChanged(obj)
            set( obj.handle , 'Enable' , CLbool2text(obj.windowmanager.isWindowOpenable(obj.windowtag) ) );
        end
        
        
        function destructor(obj)
            delete(obj.eventlistener);
            delete(obj);
        end
        
    end
    
end




