classdef CLSettingPRCWindow < CLSetting
    
    properties
        windowtype
        names
        color
        lightcolor
        derivative
        
        fighandle
        outputter
        
        locked = 0;
    end
    
    
    methods
        function obj = CLSettingPRCWindow(windowtype)
            obj = obj@CLSetting(windowtype, false, InputRestrictions.BOOL, 2, 100, 3, CLSettingsHelp.getHelp('PRC'));
            obj.windowtype = windowtype;
            
            if strcmpi(windowtype, 'prc')
                obj.names = {'PRC', 'Phase', 'Response'};
                obj.color = 'blue';
                obj.lightcolor = [0.5843    0.8157    0.9882];
                obj.derivative = 0;
                
            elseif strcmpi(windowtype, 'dprc')
                obj.names = {'dPRC', 'Phase', 'Derivative of response'};
                obj.color = 'r';
                obj.lightcolor = [0.9882 0.5843 0.5843];
                obj.derivative = 1;
                
            end
            
            obj.locked = 0;
        end
        
        function launchWindow(obj, session)
            if isempty(session)
               obj.fighandle = figure();
            else
               obj.fighandle = session.windowmanager.demandWindow(obj.names{1}); 
            end
            
            obj.outputter = GUIPRCPlot(session, obj.fighandle, obj.names, obj.color, obj.lightcolor, obj.derivative);
            set(obj.fighandle, 'Visible', 'on');
            obj.setValue(true);
        end
        
        
        function syncValue(obj)
            obj.value = ~isempty(obj.fighandle) && isvalid(obj.fighandle);
        end

        function value = getValue(obj)
            if ~obj.locked
                obj.syncValue();
                
                if ~obj.value
                    obj.fighandle = [];
                    obj.outputter = [];
                    
                end
            end
            value = obj.value;
            
            
            
        end
        
        function [valid, msg] = setValue(obj, newvalue)
            [valid, msg] = setValue@CLSetting(obj, newvalue);
            obj.locked = 0;
            
        end
        
        
       function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
           grid = cell(1, 2);
           grid{1} = GUISettingLabel(panelhandle, settings, label, options{:});
           
           grid{2} = uicontrol(panelhandle, 'style', 'pushbutton', 'String', sprintf('Open %s Plot', obj.names{1}), 'Callback', @(o,e) obj.launchWindow(session));
           box = LayoutBox(grid);
           box.widths(1) = max(suggestions.labelsize, box.widths(1)); 
           
       end
       
       
        function newobj = copy(obj, newsettings)
            newobj = CLSettingPRCWindow(obj.windowtype);
            obj.copyOver(newobj);
            newobj.locked = 1;
            newobj.value = obj.value;
        end
        
        function o = getOutputter(obj)
            if ~isempty(obj.fighandle) && isvalid(obj.fighandle) && obj.outputter.doInit()
                o = obj.outputter;
            else
                o = [];
            end
            
        end
        
    end
    
end
