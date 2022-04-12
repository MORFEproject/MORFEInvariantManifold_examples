classdef GUIBrowserPreviewPanel < handle
    
    properties
        handle
        eventlistener = {};
        layoutstructure
        
        workspace = struct();
        
        previewrenderer = [];
        
        slider = [];
    end
    
    methods
        function obj = GUIBrowserPreviewPanel(parent, model , varargin)
            obj.handle = uipanel(parent, 'Unit' , 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR ,'DeleteFcn' , @(o,e) obj.destructor(), 'ResizeFcn' , @(o,e) obj.onResize(o) ,varargin{:});
            
            obj.eventlistener{1} = model.addlistener('selectionChanged' , @(o,e) obj.selectionChanged(model));
            obj.eventlistener{2} = model.addlistener('previewTypeChanged' , @(o,e) obj.previewTypeChanged(model));

            obj.previewTypeChanged(model);
            obj.selectionChanged(model);
        end
        
        function selectionChanged(obj, previewmodel)
            obj.previewrenderer.selectionChanged(obj, previewmodel);
        end
        
        function previewTypeChanged(obj,previewmodel)
            obj.clearPanel();
            obj.previewrenderer = previewmodel.getPreviewRenderer();
            obj.previewrenderer.init(obj, previewmodel);
            
        end
        
        
        function onResize(obj,handle)
            units = get(handle, 'Units');
            set(handle,'Units' , 'Pixels');
            pos = get(handle, 'Position');
            if (~isempty(pos))
                obj.layoutstructure.makeLayoutHappen( get(handle, 'Position'));
            end
            set(handle,'Units' , units);
        end

        function clearPanel(obj)
            obj.workspace = struct();
            if ~isempty(obj.handle); delete(allchild(obj.handle)); end;
            if (~isempty(obj.layoutstructure))
                obj.layoutstructure.destructor();
                obj.layoutstructure = [];
            end
        end
        
        function destructor(obj)
        end
        
    end
    
end


function s_tableresizer(handle)
data = get(handle,'Data');
col1 = maxChar(data(:,1)) * 7 + 40;
col2 = maxChar(data(:,2)) * 7 + 40;
col3 = maxChar(data(:,3)) * 7;

set(handle , 'ColumnWidth' , {col1,col2,col3});
end
