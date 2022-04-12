classdef GUISwitchMenuItem < handle
    
    
    properties
        handle
        setter
        getter
        model
        
        eventlistener = [];
    end
    
    methods
        
        
        function obj = GUISwitchMenuItem (parent , string ,  setter , getter , model ,eventname ,  varargin )
            
            obj.handle = uimenu(parent  , 'Label' , string , 'Callback' , @(src,ev) obj.newValue()  , varargin{:} , 'DeleteFcn' , @(o,e) obj.destructor() , ...
                'Checked' , bool2str(getter()) );
            
            obj.setter = setter;
            obj.getter = getter;
            if ~isempty(model)
                obj.eventlistener = model.addlistener(eventname , @(srv,ev) obj.settingChanged());
            end
            obj.handle.UserData = obj;
        end
        
        function newValue(obj)
            
            value =  ~strcmp(get(obj.handle , 'Checked' ),'on');
            obj.setter( value );
            set(obj.handle , 'Checked' , bool2str(value));
        end
        
        function settingChanged(obj)
            set(obj.handle , 'Checked' , bool2str(obj.getter()));
        end
        
        
        function destructor(obj)
            delete(obj.eventlistener);
            delete(obj);
        end
        
    end
    
end

function result = bool2str( bool)
if  (bool)
    result = 'on';
else
    result = 'off';
end
end