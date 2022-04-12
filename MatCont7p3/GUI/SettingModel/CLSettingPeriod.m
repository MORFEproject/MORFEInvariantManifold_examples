classdef CLSettingPeriod < CLSetting
    

    methods
        function obj = CLSettingPeriod(initialvalue)
            
             obj = obj@CLSetting('Period', initialvalue, InputRestrictions.BOOL, 2, 2, Inf, CLSettingsHelp.getHelp('Period'));
        
             
        end
        
        function unselectable(obj, value)
           obj.value = value;
           obj.editable = false; 
        end
        
        function setSelectable(obj, b)
           obj.editable = b;
           if ~b
              obj.value = false; 
           end
        end
        
        function b = isSelectable(obj)
            b = obj.editable;
        end
        
        function box = renderGUI(varargin)
            box = []; %Rendering done by CLSettingParameters
        end        
        
        function newobj = copy(obj, newsettings)
           newobj = CLSettingPeriod(obj.value);
           newobj.editable = obj.editable; 
           obj.copyOver(newobj);
       end
        
    end
    
    
    
    
end