classdef GUISettingEditBox < handle
    
    properties
        handle
        eventlistener = [];
        settings
        settingname
        OKColor;
    end
    
    methods
        
         function obj = GUISettingEditBox(parent, settings, settingname, varargin )
             

             obj.settingname = settingname;
             obj.settings = settings;
            

             
            obj.handle = uicontrol(parent , 'Style' , 'edit'  , 'Unit' , 'Pixels'  , 'String' , '' , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} );  
            obj.OKColor = get(obj.handle, 'BackgroundColor');
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.setBackgroundOK();
            obj.eventlistener = settings.addlistener('settingChanged' , @(srv,ev) obj.settingChanged()); 
            obj.settingChanged();
           
        end
        
        function newValue(obj)
            string =  get(obj.handle ,'String');

            try
                if ~isempty(string)
                    x = evalin('base' , string);
                else
                    x = []; 
                end
                setting = obj.settings.getSetting(obj.settingname);
                [valid, errormsg] = setting.setValue(x);
                
                if (~valid)
                    obj.performErrorDisplay();
                    fprintf(2, sprintf('[%s] ERROR(%s): %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), obj.settingname, errormsg, string));
                    obj.settingChanged();
                else
                    obj.settings.refresh();
                end
                
            catch error
          
                fprintf(2, sprintf('[%s] ERROR(%s): %s: %s, value: %s\n\n',datetime('now', 'format', 'HH:mm:ss'), obj.settingname, error.message, string));
                obj.performErrorDisplay();
                obj.settingChanged(); %restore original value
            end
       
  
        end
        function performErrorDisplay(obj)
            obj.setBackgroundERROR();
            pause(0.3)
            obj.setBackgroundOK();
            
        end
        function setBackgroundOK(obj)
           set(obj.handle, 'BackgroundColor' , [1 1 1]);
           set(obj.handle, 'BackgroundColor' , obj.OKColor);
        end

        function setBackgroundERROR(obj)
           set(obj.handle, 'BackgroundColor' ,  [1    0.3    0.3]);
        end
                       
        
        function settingChanged(obj)
           set(obj.handle , 'String' , obj.settings.getSetting(obj.settingname).toString()); 
        end
        
        function destructor(obj)
            delete(obj.eventlistener);
            delete(obj);
        end
        
        function e = Extent(obj)
           e = obj.handle.Extent; 
        end
        function set(obj, varargin)
           %disp(varargin{2});
           set(obj.handle, varargin{:}); 
        end
    end
    
end

