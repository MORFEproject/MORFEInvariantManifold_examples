classdef CLSettingCustomParameters < CLSetting
    
    properties
        parameters
        parametersfree
        parametersVisible
        selectFree
    end
    
    methods
        function obj = CLSettingCustomParameters(settings, parameters, selectFree)
            
            dim = length(parameters);
            obj = obj@CLSetting('customparam', zeros(1, dim), InputRestrictions.vector(dim), 2, 3, 100, '~');
            
            obj.parameters = parameters;
            obj.parametersfree = zeros(1, dim);
            obj.parametersVisible = 1;
            
            if ~isempty(settings)
                obj.installGhosts(settings, selectFree);
            end
            obj.selectFree = selectFree;
        end
        
        function list = getActive(obj)
            list = find(obj.parametersfree);
        end

        
        
        function revive(obj, settings, selectFree, visible)
            if nargin < 4
               visible = true; 
            end

            obj.selectFree = selectFree;
            obj.setVisible(true);
            if ~isempty(settings)
                obj.installGhosts(settings, selectFree);
            end
            obj.parametersVisible = visible;
        end
        
        function installGhosts(obj, settings, selectFree)
            for i = 1:length(obj.parameters)
                parameter = obj.parameters{i};
                settings.addSetting(['custom_', parameter], CLSettingParameter(obj, i, InputRestrictions.POS));
                if selectFree; settings.addSetting(['custom_', parameter, '_select'], CLSettingParameterSelect(obj, i, 'select', length(obj.parameters))); end
            end
        end
        
        function newobj = copy(obj, newsettings)
            newobj = CLSettingCustomParameters(newsettings, obj.parameters, obj.selectFree);
            newobj.parametersfree = obj.parametersfree;
            newobj.parametersVisible = obj.parametersVisible;
            newobj.value = obj.value;
            obj.copyOver(newobj);
        end
        
        
        function val = get(obj, tag, index)
            val = obj.parametersfree(index);
        end
        
        function set(obj,tag, index, newval)
            obj.parametersfree(index) = newval;
        end
        
        
        
        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            
            selectiongrid = cell(length(obj.parameters), obj.selectFree + 1);
            grid = cell(length(obj.parameters), 2);
            
            for i = 1:length(obj.parameters)
                    parameter = obj.parameters{i};
                    index = 1;
                    if obj.selectFree
                        selectiongrid{i, index} = GUISettingSwitch('radiobutton', panelhandle, ...
                            settings, ['custom_', parameter, '_select'], [options, {'TooltipString'}, {sprintf('select "%s" as active parameter', parameter)}]);
                        index = index + 1;
                    end
                    selectiongrid{i, index} = uicontrol( panelhandle, 'Style', 'text', 'String', parameter,'HorizontalAlignment', 'left' , options{:}, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                    
                    subbox = LayoutBox(selectiongrid(i, :));
                    
                    grid{i, 1} = subbox;
                    if obj.parametersVisible
                        grid{i, 2} = GUISettingEditBox(panelhandle, settings, ['custom_', parameter], options{:}, 'Enable', CLbool2text(obj.editable));
                    else
                        grid{i, 2} = [];
                    end
            end
            
           
            box1 = LayoutBox(grid);
            box1.widths(1) = max(suggestions.labelsize, box1.widths(1)); 
 
            if ~isempty(session)
                curvelabel = session.getCompConf().getCurveLabel();
                if contains(curvelabel, 'Het')
                   header = 'Heteroclinic parameters';
                else
                    header = 'Homoclinic parameters';
                end
            else
               header = 'Connecting Orbit parameters'; 
            end
            
            sectionsettings = DefaultValues.SECTIONNAMESETTING;
            grid = {GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10); uicontrol(panelhandle, 'Style' , 'Text', 'String', header, sectionsettings{:}); box1};
            box = LayoutBox(grid);
        end
        
        
        
    end
end
