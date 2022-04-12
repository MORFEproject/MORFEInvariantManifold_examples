classdef GUIPreviewRendererSystems < GUIPreviewRenderer
    
    methods
        
        function init(~, previewpanel, ~)
            
            mainnode = LayoutNode(-1,-1,'vertical');
            previewpanel.layoutstructure = mainnode;
            position = get(previewpanel.handle , 'Position');
            
            previewpanel.workspace.tablehandle = uitable(previewpanel.handle ,'Unit', 'Pixels');
            mainnode.addHandle( 1, 1 , uicontrol(previewpanel.handle, 'style', 'text', 'String', ''));
            mainnode.addHandle( 99, 99 , previewpanel.workspace.tablehandle , 'minsize' , [Inf,Inf]);
            previewpanel.layoutstructure.makeLayoutHappen(position);
            
            
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            system = previewmodel.getInfoObject();
            
            if isempty(system); return; end
            
            tablehandle = previewpanel.workspace.tablehandle;

            coords = system.getCoordinates();
            params = system.getParameters();
            derivatives = system.getDerInfo();
            equations = system.getEquations();
            
            table = cell( 1+1+1+1+1+ size(equations,1) , 2);
            table{1,1}  = 'Name';
            table{2,1} = 'Time';
            table{3,1} = 'Coordinates';
            table{4,1} = 'Parameters';
            table{5,1} = 'Derivatives';
            table{6,1} = 'Userfunctions';
            table{7,1} = 'Equations';

            table{1,2} = system.getName();
            table{2,2} = system.getTimeName();
            table{3,2} = cellarray2singlestr(coords);
            table{4,2} = cellarray2singlestr(params);
            table{5,2} = derivatives;
            table{6,2} = system.getUFString();
            for i = 1:size(equations,1)
                table{6+i , 2} = equations(i,:);
            end
            col1w = maxChar(table(:,1)) * DefaultValues.LETTERDIMENSION(1);
            col2w = maxChar(table(:,2)) * DefaultValues.LETTERDIMENSION(1);

            set(tablehandle, 'Units' , 'pixels' , 'ColumnWidth' , {col1w , col2w} ,  'Data' , table , 'RowName', [] , 'ColumnName' , [] );            
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
function str = cellarray2singlestr(array)
str = '';
if (~isempty(array))
    str = array{1};
end
for i = 2:length(array)
    str = [str ', ' array{i}];
end
end


