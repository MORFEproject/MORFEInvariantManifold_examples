classdef GUISettingPopupmenu < handle
    
    properties
        handle
        eventlistener = [];
        settings
        settingname
    end
    
    methods
        
        function obj = GUISettingPopupmenu(parent, settings, settingname, varargin )
            
            
            obj.settingname = settingname;
            obj.settings = settings;
            setting = obj.settings.getSetting(settingname);
            strings = setting.validitycheck.categories;
            
            obj.handle = uicontrol(parent , 'Style' , 'popupmenu'  , 'Unit' , 'Pixels'  , 'String' , strings, 'Value', 1 , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} );
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.eventlistener = settings.addlistener('settingChanged' , @(srv,ev) obj.settingChanged());
            obj.settingChanged();
            
        end
        
        function newValue(obj)
            strings =  get(obj.handle ,'String');
            value = strings{get(obj.handle, 'Value')};
            obj.settings.setValue(obj.settingname, value);
   
        end
        


       
        
        function settingChanged(obj)
          value =   obj.settings.getValue(obj.settingname);
          strings =  get(obj.handle ,'String');
          set(obj.handle, 'Value', find(strcmp(strings, value)))
            
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

