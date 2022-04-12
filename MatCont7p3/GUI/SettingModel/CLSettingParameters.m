classdef CLSettingParameters < CLSetting
    
    properties
        parameters
        parametersfree
        parametersbranch
        selectFree
        selectBranch
        
        editFree = true;
        editBranch = true;
    end
    
    methods
        function obj = CLSettingParameters(settings, parameters, selectFree, selectBranch)
            dim = length(parameters);
            obj = obj@CLSetting('param', zeros(1, dim), InputRestrictions.vector(dim), 2, 2, 0, '~~~');
            obj.parameters = parameters;
            obj.parametersfree = zeros(1, dim);
            obj.parametersbranch = zeros(1, dim);
            
            if nargin < 3
                selectFree = false;
            end
            if nargin < 4
               selectBranch = false; 
            end
            if ~isempty(settings)
                obj.installGhosts(settings, selectFree, selectBranch);
            end
            obj.selectFree = selectFree;
            obj.selectBranch = selectBranch;
        end
        function list = getActive(obj)
            list = find(obj.parametersfree);
        end
        
        function b =  hasActive(obj)
            if ~obj.selectFree
                b = [];
            else
               b = obj.getActive(); 
            end
            
        end
        
        
        function list = getBranch(obj)
            if obj.selectBranch
                list = find(obj.parametersbranch);
            else
                list = [];
            end
        end
        
        function setEditable(obj, edit, editFree)
           obj.editFree = editFree;
           obj.editable = edit;
        end
        
        function revive(obj,settings, selectFree, selectBranch)
            if nargin < 4
               selectBranch = false; 
            end
            obj.selectFree = selectFree;
            obj.selectBranch = selectBranch;
            obj.setVisible(true);
            obj.editable = true;
            obj.editFree = true;
            if ~isempty(settings)
                obj.installGhosts(settings, selectFree, selectBranch);
            end
            
        end
        
        function installGhosts(obj, settings, selectFree, selectBranch)
            for i = 1:length(obj.parameters)
                parameter = obj.parameters{i};
                settings.addSetting(['pa_', parameter], CLSettingParameter(obj, i));
                if selectFree; settings.addSetting(['pa_', parameter, '_select'], CLSettingParameterSelect(obj, i, 'select', length(obj.parameters))); end
                if selectBranch; settings.addSetting(['pa_', parameter, '_branch'], CLSettingParameterSelect(obj, i, 'branch', 2*length(obj.parameters))); end
            end
        end
        
        function newobj = copy(obj, newsettings)
            newobj = CLSettingParameters(newsettings, obj.parameters, obj.selectFree, obj.selectBranch);
            newobj.parametersfree = obj.parametersfree;
            newobj.parametersbranch = obj.parametersbranch;
            newobj.value = obj.value;
            obj.copyOver(newobj);
        end
        
        
        function val = get(obj, tag, index)
            if strcmp(tag, 'select')
                val = obj.parametersfree(index);
            else
                val = obj.parametersbranch(index);
            end
        end
        function set(obj,tag, index, newval)
            if strcmp(tag, 'select')
                obj.parametersfree(index) = newval;
            else
                obj.parametersbranch(index) = newval;
            end
        end
        
        
        
        function box = renderGUI(obj, session, settings, label, panelhandle, options, suggestions)
            selectiongrid = cell(length(obj.parameters), obj.selectFree+obj.selectBranch+1);
            grid = cell(length(obj.parameters), 2);
            len = obj.selectFree+obj.selectBranch;
            
            for i = 1:length(obj.parameters)
                parameter = obj.parameters{i};
                index = 1;
                if obj.selectBranch
                    selectiongrid{i, index} = GUISettingSwitch('radiobutton', panelhandle, ...
                        settings, ['pa_', parameter, '_branch'],  [options, {'Enable'}, {CLbool2text(obj.editBranch)}, {'TooltipString'}, {sprintf('select "%s" as branch parameter', parameter)}]);
                    index = index + 1;
                end
                if obj.selectFree
                    selectiongrid{i, index} = GUISettingSwitch('radiobutton', panelhandle, ...
                        settings, ['pa_', parameter, '_select'], [options, {'Enable'}, {CLbool2text(obj.editFree)}, {'TooltipString'}, {sprintf('select "%s" as active parameter', parameter)}]);
                    index = index + 1;
                end
                selectiongrid{i, index} = uicontrol( panelhandle, 'Style', 'text', 'String', parameter,'HorizontalAlignment', 'left' , options{:}, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                
                subbox = LayoutBox(selectiongrid(i, :));
                
                grid{i, 1} = subbox;
                grid{i, 2} = GUISettingEditBox(panelhandle, settings, ['pa_', parameter], options{:}, 'Enable', CLbool2text(obj.editable));
                
            end
            
            period = settings.getSetting('Period');
            if ~isempty(period) && period.isVisible()
                selectiongrid = cell(1, obj.selectFree+obj.selectBranch+1);
                index = 1;
                if obj.selectBranch
                    %selectiongrid{index} = uicontrol('style', 'radiobutton', 'Visible', 'off');
                    selectiongrid{index} = [];
                    index + 1;
                end
                
                if obj.selectFree
                    selectiongrid{index} = GUISettingSwitch('radiobutton', panelhandle, settings, 'Period', [options, {'Visible'}, {CLbool2text(period.isSelectable())}]);
                    index = index + 1;                    
                end    
                
                selectiongrid{index} = uicontrol( panelhandle, 'Style', 'text', 'String', 'Period','HorizontalAlignment', 'left' , options{:}, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                
                subbox = LayoutBox(selectiongrid);
                
                i = length(obj.parameters) + 1;
                grid{i, 1} = subbox;
                
                ip = settings.getSetting('IP');
                if isfield(ip.currentpointdata, 'data') && isfield(ip.currentpointdata.data, 'T')
                    value = num2str( ip.currentpointdata.data.T, '%.14g');
                    visible = 'on';
                else
                    visible = 'off'; value = '~';
                end
                grid{i, 2} = uicontrol(panelhandle, 'Style', 'Text', 'String', value, options{:}, 'Visible', visible, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);                
            end
            box = LayoutBox(grid);
            box.widths(1) = max(suggestions.labelsize, box.widths(1)); 
            %GUIrender = obj.validitycheck.getGUIRenderer();
            %box2 = LayoutBox({GUIrender(panelhandle, settings, label, options{:} , 'BackgroundColor', [0.90 0.90 0.90] , 'Enable', CLbool2text(obj.editable))  });
            %box = LayoutBox({box1; box2});

        end
        
        
        function b = sanityCheck(obj, settings)
            b = obj.internalSanityCheck(settings);
            
            if ~b
                for i = 1:length(obj.parameters)
                    parameter = obj.parameters{i};
                    settings.removeSetting(['pa_', parameter]);
                    settings.removeSetting(['pa_', parameter '_select']);
                    settings.removeSetting(['pa_', parameter '_branch']);
                end
            end
            
        end
        
        function b = internalSanityCheck(obj, settings)
            system = settings.system;
            if isempty(system); b = 0; return; end
            parameters = system.getParameters();
            if length(parameters) ~= length(obj.parameters); b = 0; return; end
            
            for k = 1:length(obj.parameters)
                if ~strcmp(obj.parameters{k}, parameters{k})
                    b = 0;
                    return
                end
            end
            b = 1;
        end
        
    end
end
