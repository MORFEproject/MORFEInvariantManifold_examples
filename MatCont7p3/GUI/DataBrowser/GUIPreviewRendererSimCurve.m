classdef GUIPreviewRendererSimCurve < GUIPreviewRenderer
    
    methods
        
        function init(~, previewpanel, previewmodel)
            
            grid = cell(2*3, 1);
            solution =  previewmodel.getPreviewData();
            system = solution.settings.system;           
            
            
            
            
            previewpanel.workspace.npointhandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', [],'Data', num2cell(pi*ones(1, 1)) );
            grid{1} = previewpanel.workspace.npointhandle;
                        
            coordinates = [system.getTimeName(), system.getCoordinates()];
            
            grid{2} = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'String', 'coordinates', 'HorizontalAlignment', 'left', 'BackgroundColor', DefaultValues.BACKGROUNDCOLOR);
            
            previewpanel.workspace.coordtablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', coordinates, 'Data', num2cell(pi*ones(length(coordinates), 1)));
            grid{3} = previewpanel.workspace.coordtablehandle;
            
            
            parameters = system.getParameters();
            grid{4} = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'String', 'parameters', 'HorizontalAlignment', 'left', 'BackgroundColor', DefaultValues.BACKGROUNDCOLOR);
            
            previewpanel.workspace.paramtablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels', 'ColumnName', [], 'RowName', parameters, 'Data', num2cell(pi*ones(length(parameters), 1)));
            grid{5} = previewpanel.workspace.paramtablehandle;
            
            previewpanel.workspace.msglabel = uicontrol(previewpanel.handle, 'Style', 'text', 'Unit', 'Pixels', 'horizontalalignment', 'left', 'String' , ['this' newline 'is a warning'], 'BackgroundColor', DefaultValues.BACKGROUNDCOLOR);
            grid{6} = previewpanel.workspace.msglabel;
            
            box = LayoutBox(grid);
            box.sourcepanel = previewpanel;
            
            box.makeLayoutHappen(get(previewpanel.handle , 'Position'));
            previewpanel.layoutstructure = box;           
            
            
            
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            
            data = previewmodel.getInfoObject();
            if isfield(data, 'npoint')
                set(previewpanel.workspace.npointhandle, 'RowName', {'npoint'}, 'Data', {sprintf('%.8g',data.npoint)});
                set(previewpanel.workspace.msglabel, 'String', '');
            else
                set(previewpanel.workspace.npointhandle, 'RowName', {'interval'}, 'Data', {sprintf('%.8g',data.interval)});
                set(previewpanel.workspace.msglabel, 'String', data.selectmsg);
            end
            set(previewpanel.workspace.coordtablehandle, 'Data', arrayfun(@(x) sprintf('%.8g',x), [data.time; data.coord(:)], 'un', 0));
            set(previewpanel.workspace.paramtablehandle, 'Data', arrayfun(@(x) sprintf('%.8g',x), data.param(:), 'un', 0));
            tableresizer(previewpanel.workspace.coordtablehandle);
            tableresizer(previewpanel.workspace.paramtablehandle);
            
        end
    end
    

    methods(Static)

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
col1 = maxChar(data(:,1)) * 10 + 40;
set(handle , 'ColumnWidth' , {col1});
end

