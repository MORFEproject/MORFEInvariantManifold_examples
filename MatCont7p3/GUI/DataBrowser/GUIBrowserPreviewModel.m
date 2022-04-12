classdef GUIBrowserPreviewModel < handle
    
    properties
        %session
        
        browsermodel
        current
        
        eventlistener = {}
        info = [];
        
        index = -1;
    end
    
    events
        previewTypeChanged
        selectionChanged
    end
    
    methods
        function obj = GUIBrowserPreviewModel(session, browsermodel)
            %obj.session = session;
            obj.browsermodel = browsermodel;
            obj.eventlistener{1} = browsermodel.addlistener('indexChanged', @(o,e) obj.onIndexChanged(obj.current));
            obj.eventlistener{2} = browsermodel.addlistener('listChanged', @(o,e) obj.onListChanged(browsermodel));
            
            obj.onListChanged(browsermodel);
        end
        
        function onIndexChanged(obj,currentmodel)
            obj.info = currentmodel.getSelectedPreviewData();
            if currentmodel.isValidSelection()
                obj.notify('selectionChanged');
            end
        end

        function onListChanged(obj,browsermodel)
            obj.current = browsermodel.getCurrent();
            obj.notify('previewTypeChanged');
            obj.onIndexChanged(obj.current);
        end    
        
        function info = getInfoObject(obj)
            info = obj.info;
        end
        
        
        function renderer = getPreviewRenderer(obj)
            renderer = obj.current.getPreviewRenderer();
        end
        function setSelectedIndex(obj, index)
            obj.index = index;
        end
        
        function i = getSelectedIndex(obj)
            i = obj.index;
        end
        
        function list = getSelectList(obj)
            list = obj.info;
        end
        function n = getNrElements(obj)
            n = length(obj.info);
        end
        function unSetIndex(obj)
            obj.index = -1;
        end
        function selectItem(obj, index)
        end
        function goUp(obj)
        end
        function it = getSelectedItemName(obj)
            it = '';
        end
        function o = getPreviewData(obj)
           o = obj.browsermodel.getCurrent().getPreviewData(); 
        end
    end
    
end

