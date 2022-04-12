classdef GUIPreviewRendererContCurve < GUIPreviewRenderer
    
    methods
        
        function init(~, previewpanel, previewmodel)
            grid = cell(2*3+1, 1);
            solution =  previewmodel.getPreviewData();
            system = solution.settings.system;
            settings = solution.settings;
            
            previewpanel.workspace.npointhandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', [],'Data', num2cell(pi*ones(1, 1)) );
            grid{1} = previewpanel.workspace.npointhandle;
            
            if settings.getVal('ntst', 0) == 0 %not a cycle.
                coordinates = system.getCoordinates();
                grid{2} = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'String', 'coordinates', 'HorizontalAlignment', 'left');
                previewpanel.workspace.coordtablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', coordinates,  'Data', num2cell(pi*ones(length(coordinates), 1)));
                grid{3} = previewpanel.workspace.coordtablehandle;
                previewpanel.workspace.period = false;
            else
                grid{2} = [];
                previewpanel.workspace.coordtablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', {'Period'},  'Data', num2cell(pi*ones(1, 1)));
                grid{3} = previewpanel.workspace.coordtablehandle;
                previewpanel.workspace.period = true;
                
            end
            parameters = system.getParameters();
            grid{4} = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'String', 'parameters', 'HorizontalAlignment', 'left');
            previewpanel.workspace.paramtablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', parameters,  'Data', num2cell(pi*ones(length(parameters), 1)));
            grid{5} = previewpanel.workspace.paramtablehandle;
            
            previewpanel.workspace.bifdatalabel = uicontrol(previewpanel.handle,'Style' , 'text', 'Visible', 'on', 'Unit', 'Pixels', 'String', 'non-degeneracy conditions', 'HorizontalAlignment', 'left');
            previewpanel.workspace.biftablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'Visible', 'off', 'ColumnName', [], 'RowName', {'tmp', 'tmp2'}, 'Data', num2cell(pi*ones(2, 1)));
            grid{6} = previewpanel.workspace.bifdatalabel;
            grid{7} = previewpanel.workspace.biftablehandle;
            

            box = LayoutBox(grid);
            box.sourcepanel = previewpanel;
            
            box.makeLayoutHappen(get(previewpanel.handle , 'Position'));
            previewpanel.layoutstructure = box;
            
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            data = previewmodel.getInfoObject();
            
            set(previewpanel.workspace.npointhandle, 'RowName', {'npoint'}, 'Data', {sprintf('%.8g',data.npoint)});
            
            if previewpanel.workspace.period
                set(previewpanel.workspace.coordtablehandle, 'Data', {sprintf('%.8g',data.period)});
            else
                set(previewpanel.workspace.coordtablehandle, 'Data', arrayfun(@(x) sprintf('%.8g',x), data.coord(:), 'un', 0));
            end
            set(previewpanel.workspace.paramtablehandle, 'Data', arrayfun(@(x) sprintf('%.8g',x), data.param(:), 'un', 0));
            
            tableresizer(previewpanel.workspace.coordtablehandle);
            tableresizer(previewpanel.workspace.paramtablehandle);
            
            tablecontent = GUIPreviewRendererContCurve.getBifurcationData(data.detect, @(x) strip(strip(sprintf('%.8e, ',x)), ','));
            %{
            if isempty(tablecontent)
                set(previewpanel.workspace.bifdatalabel, 'Visible', 'off');
                set(previewpanel.workspace.biftablehandle, 'visible', 'off');
                
            else
            %}
            %set(previewpanel.workspace.bifdatalabel, 'Visible', 'on');
            if isempty(tablecontent)
                set(previewpanel.workspace.biftablehandle, 'visible', 'on', 'RowName', [], 'Data', []);
            else
                set(previewpanel.workspace.biftablehandle, 'visible', 'on', 'RowName', tablecontent(:, 1), 'Data', tablecontent(:, 2), 'ColumnWidth', {Inf});
            end
            tableresizer(previewpanel.workspace.biftablehandle);
            %end
            
            
        end
    end
    
    
    methods(Static)
        function content = getBifurcationData(s, tostring)
            content = {};
            if isempty(s) || ~isfield(s, 'data'); return; end
            pointtype = strip(s.label);
            if isfield(s, 'internallabel') && ~isempty(s.internallabel)
                pointtype = strip(s.internallabel);
            end
            
            if strcmp(pointtype, 'H') && isfield(s.data, 'lyapunov')
                content = {'1st Lyapunov', tostring(s.data.lyapunov)};
            elseif  strcmp(pointtype, 'LP') && isfield(s.data, 'a')
                content = {'a', tostring(s.data.a)};
            elseif strcmp(pointtype, 'LPC') && isfield(s.data, 'lpccoefficient')
                content = {'lpccoefficient', tostring(s.data.lpccoefficient)};
            elseif strcmp(pointtype, 'NS') && isfield(s.data, 'nscoefficient')
                content = {'nscoefficient', tostring(s.data.nscoefficient)};
            elseif strcmp(pointtype, 'PD') && isfield(s.data, 'pdcoefficient')
                content = {'pdcoefficient', tostring(s.data.pdcoefficient)};
            elseif strcmp(pointtype, 'CPC') && isfield(s.data, 'cpccoefficient')
                content = {'cpccoefficient', tostring(s.data.cpccoefficient)};
            elseif strcmp(pointtype, 'LPPD') && isfield(s.data, 'ffcoefficients')
                content = {'ffcoefficients', tostring(s.data.ffcoefficients)};
            elseif strcmp(pointtype, 'GPD') && isfield(s.data, 'gpdcoefficient')
                content = {'gpdcoefficient', tostring(s.data.gpdcoefficient)};
            elseif isfield(s.data, 'c')
                content = {'c', tostring(s.data.c)};
            end
            
        end
    end
       
    
    
end





function m = maxChar(cellstr)
l = zeros(1, length(cellstr));
for i = 1:length(cellstr)
    l(i) = length(cellstr{i});
end
m = max(l);
end


function tableresizer(handle)
data = get(handle,'Data');
if ~isempty(data)
col1 = maxChar(data(:,1)) * 10 + 40;
set(handle , 'ColumnWidth' , {col1});
end
end

