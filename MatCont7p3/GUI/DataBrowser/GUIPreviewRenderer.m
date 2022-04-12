classdef GUIPreviewRenderer
    
    methods
    
        function init(~, previewpanel, previewmodel)
            previewpanel.layoutstructure = LayoutNode(-1, -1, 'vertical');
            fprintf(2, 'rendering init\n');
        end
        
        function selectionChanged(~, previewpanel, previewmodel)
            fprintf(2, 'rendering selection changed.\n');
        end
    
    end
    
    
end
