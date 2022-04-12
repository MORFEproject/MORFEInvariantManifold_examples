classdef GUIPreviewRendererSolutions < GUIPreviewRenderer
    
    methods
        
        function init(~, previewpanel, ~)
            mainnode = LayoutNode(-1,-1,'vertical');
            mainnode.setOptions('Add',true);
            previewpanel.layoutstructure = mainnode;
            previewpanel.workspace.typename = '~';
            
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            solution =  previewmodel.getInfoObject();
            
            if ~isempty(solution)
                renderer = solution.compbranch.getOutputInterpreter().getPreviewPanelRenderer();
                renderer(solution, previewpanel, previewmodel);
            end
            
        end
    end
    methods(Static)
        
        
        function renderContinuationCurve(solution, previewpanel, previewmodel)
            
            if ~strcmp(previewpanel.workspace.typename, 'continuationcurve')
                set(previewpanel.handle, 'Visible', 'off');
                %fprintf(2, 'reconfigure Solutionpreviewpanel\n');
                previewpanel.workspace.typename = 'continuationcurve';
                
                previewpanel.layoutstructure.makeEmpty();
                mainnode = previewpanel.layoutstructure;
                
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, ...
                    'String' , 'Npoints' , 'HorizontalAlignment' , 'left'  ),'halign' , 'l');
                previewpanel.workspace.npointshandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  'String' , '                         ');
                subnode.addHandle(1,2,previewpanel.workspace.npointshandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  ...
                    'String' , 'Initial Pointtype', 'HorizontalAlignment' , 'left'  ),'halign' , 'l');
                previewpanel.workspace.pthandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'String' , '                         ');
                subnode.addHandle(1,2,previewpanel.workspace.pthandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, ...
                    'String' , 'Curvetype' , 'HorizontalAlignment' , 'left' ),'halign' , 'l');
                previewpanel.workspace.cthandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  'String' , '                          ');
                subnode.addHandle(1,2,previewpanel.workspace.cthandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                
                previewpanel.workspace.table = uitable(previewpanel.handle,'Unit', 'Pixels', 'RowName' , [] , 'ColumnName' , {'Type' , 'Index' , 'Message'});
                mainnode.addHandle(7 , 1 , previewpanel.workspace.table , 'minsize' , [Inf,Inf]);
                
                
                previewpanel.layoutstructure.makeLayoutHappen(  get(previewpanel.handle , 'Position') );
                set(previewpanel.handle, 'Visible', 'on');
            end
            set(previewpanel.workspace.npointshandle , 'String' , num2str(size(solution.x,2)));
            set(previewpanel.workspace.pthandle , 'String' , solution.settings.getValue('IP').getLabel());
            set(previewpanel.workspace.cthandle , 'String' , solution.compbranch.toString());
            set(previewpanel.workspace.table , 'Data' , sToCell(solution.s));
            s_tableresizer(previewpanel.workspace.table);
            
        end
        function renderSimulation(solution, previewpanel, previewmodel)
            if ~strcmp(previewpanel.workspace.typename, 'simulationcurve')
                previewpanel.workspace.typename = 'simulationcurve';
                set(previewpanel.handle, 'Visible', 'off');
                
                previewpanel.layoutstructure.makeEmpty();
                mainnode = previewpanel.layoutstructure;
                
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, ...
                    'String' , 'Npoints' , 'HorizontalAlignment' , 'left'  ),'halign' , 'l');
                previewpanel.workspace.npointshandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  'String' , '                         ');
                subnode.addHandle(1,2,previewpanel.workspace.npointshandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  ...
                    'String' , 'Initial Pointtype', 'HorizontalAlignment' , 'left'  ),'halign' , 'l');
                previewpanel.workspace.pthandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'String' , '                         ');
                subnode.addHandle(1,2,previewpanel.workspace.pthandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, ...
                    'String' , 'Integrator' , 'HorizontalAlignment' , 'left' ),'halign' , 'l');
                previewpanel.workspace.cthandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  'String' , '                          ');
                subnode.addHandle(1,2,previewpanel.workspace.cthandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                subnode = LayoutNode(1,1);
                subnode.addHandle(1,1, uicontrol(previewpanel.handle,'Style' , 'text', 'Unit', 'Pixels', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, ...
                    'String' , 'tspan' , 'HorizontalAlignment' , 'left' ),'halign' , 'l');
                previewpanel.workspace.tspanhandle = uicontrol(previewpanel.handle,'Style' , 'text', 'Unit','Pixels','BackgroundColor' , DefaultValues.BACKGROUNDCOLOR,  'String' , '                          ');
                subnode.addHandle(1,2,previewpanel.workspace.tspanhandle,'halign' , 'l', 'minsize', [Inf, DefaultValues.LETTERDIMENSION(2)]);
                mainnode.addNode(subnode);
                %
                
                
                previewpanel.layoutstructure.makeLayoutHappen(  get(previewpanel.handle , 'Position') );
                set(previewpanel.handle, 'Visible', 'on');
            end
            set(previewpanel.workspace.npointshandle , 'String' , num2str(size(solution.y,1)));
            set(previewpanel.workspace.pthandle , 'String' , solution.settings.getValue('IP').getLabel());
            set(previewpanel.workspace.cthandle , 'String' , solution.compbranch.toString());
            set(previewpanel.workspace.tspanhandle , 'String' , sprintf('[%.5g, %.5g]',solution.t(1), solution.t(end)) );
            
            
            
        end
    end
    
    
    
    
end

function c = sToCell(s)
c = cell(length(s), 3);
for i = 1:length(s)
    c{i, 1} = strip(s(i).label);
    c{i, 2} = s(i).index;
    c{i, 3} = s(i).msg;
end


end

function m = maxChar(cellstr)
l = zeros(1, length(cellstr));
for i = 1:length(cellstr)
    l(i) = length(cellstr{i});
end
m = max(l);
end


function s_tableresizer(handle)
data = get(handle,'Data');
col1 = maxChar(data(:,1)) * 10 + 40;
col2 = maxChar(data(:,2)) * 10 + 40;
col3 = maxChar(data(:,3)) * 10;

set(handle , 'ColumnWidth' , {col1,col2,col3});
end

