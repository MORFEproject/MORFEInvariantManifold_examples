classdef GUISettingSwitch < handle
    
    properties
        handle
        eventlistener = [];
        settings
        settingname
        
        extentpos
    end
    
    methods
        
        function obj = GUISettingSwitch(style, parent, settings, settingname,  varargin )
            obj.settingname = settingname;
            obj.settings = settings;
            obj.handle = uicontrol(parent , 'Style' , style, 'Unit' , 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , 'Value' , 0 , 'Callback' , @(src,ev) obj.newValue()  , varargin{1}{:} );
            
            obj.extentpos = obj.handle.Extent;
            obj.extentpos(3) = 4*obj.extentpos(3);
            
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            obj.eventlistener = settings.addlistener('settingChanged' , @(srv,ev) obj.settingChanged());
            obj.settingChanged();
            
        end
        
        function newValue(obj)
            b = obj.settings.setValue(obj.settingname, get(obj.handle, 'Value'));
        end
        
        function settingChanged(obj)
            set(obj.handle , 'Value' , obj.settings.getSetting(obj.settingname).getValue());
        end
        
        function destructor(obj)
            delete(obj.eventlistener);
            delete(obj);
        end
        
        function e = Extent(obj)
            e = obj.extentpos;
        end
        function set(obj, varargin)
            set(obj.handle, varargin{:});
        end
    end
    
end

