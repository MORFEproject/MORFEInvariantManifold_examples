classdef CLSettingEventFunction < CLSettingInterface
    
    
    properties
        functionhandle; 
        visible;
        displayname;
    end
    
    
    methods
        function obj = CLSettingEventFunction(name)
            obj.functionhandle = [];
            obj.visible = 1;
            obj.displayname = name;
        end
        
        function value = getValue(obj)
            value = obj.functionhandle;
        end
        
        function [valid, msg] = setValue(obj, newvalue)
             valid = 1; msg = '';
             
            
            try
                value = strip(newvalue, '@');
                if isempty(value)
                    obj.functionhandle = [];
                else
                    obj.functionhandle = evalin('base', ['@' value]);
                end
                
            catch error
                valid = 0;
                msg = error.message;
                obj.functionhandle = [];
            end
        end
        function id = getGroupID(~); id = 3; end
        function id = getSubGroupID(~); id = 1; end
        function id = getItemID(~); id = 2; end        
        
        function s = toString(obj)
            if isempty(obj.functionhandle)
                s = '<disabled>';
            else
                s = func2str(obj.functionhandle);
            end
        end
        
        function  b = isVisible(obj)
           b = obj.visible;
        end
        function setVisible(obj, b)
           obj.visible = b; 
        end
        
        function h = getHelpStr(~); h = '~ ~ ~'; end
        
        function t = getValueType(~); t = 'STRING'; end
        
        function newobj = copy(obj, ~)
            newobj = CLSettingEventFunction(obj.displayname);
            newobj.visible = obj.visible;
            newobj.functionhandle = obj.functionhandle;
        end       

        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            grid = cell(1, 2);
            grid{1} = GUISettingLabel(panelhandle, settings, label, options{:});
            grid{2} = uicontrol(panelhandle', 'style', 'edit',  'String', obj.toString(), 'Callback' , @(src,ev) obj.callbackChange(src), options{:});
            box = LayoutBox(grid);
            box.widths(1) = max(suggestions.labelsize, box.widths(1));
        end
        
        function callbackChange(obj, handle)
            str = get(handle, 'String');
            
            
            [valid, msg] = obj.setValue(str);
            if ~valid
                fprintf(2, 'error: %s\n\n', msg);
            end
            
            set(handle, 'String', obj.toString());
            
        end
    end
end




