classdef GUIPreviewRendererDiagrams < GUIPreviewRenderer
    
    methods
        
        function init(~, previewpanel, ~)
            mainnode = LayoutNode(-1,-1,'vertical');
            previewpanel.layoutstructure = mainnode;
            previewpanel.workspace.listhandle = uicontrol(previewpanel.handle,'Unit', 'Pixels', 'Style', 'edit' , 'BackgroundColor' , 'white' , 'Max' , 10 , 'Enable' , 'inactive' );
            mainnode.addHandle( 1, 1 , previewpanel.workspace.listhandle , 'minsize' , [Inf,Inf]);
            previewpanel.layoutstructure.makeLayoutHappen(  get(previewpanel.handle , 'Position'));
            
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            diagram = previewmodel.getInfoObject();
            set(previewpanel.workspace.listhandle, 'String' , diagram.getCleanSolutionNames());
        end
        
    end
    
    
end
