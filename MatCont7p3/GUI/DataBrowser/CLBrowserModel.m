classdef CLBrowserModel < handle

    properties
        stack
        current
        
        %session
        
    end
    
    events
        listChanged
        indexChanged
        shutdown
    end
    
    methods
        function obj = CLBrowserModel()
           % obj.session = session;
            obj.stack = {};
            obj.current = []; 
        end
        
        function initCurrentModel(obj, model)
           obj.stack = {};
           obj.current = model;
           obj.current.initIndex();
        end
        
        
        function list = getSelectList(obj)
            list = obj.current.getList();
        end
        
        function n = getNrElements(obj)
            n = length(obj.current.getList());
        end
        
        
        function setSelectedIndex(obj, index)
            obj.current.index = index;
            obj.notify('indexChanged');
        end
        
        function fixIndex(obj)
           if isempty(obj.getSelectList())
               obj.current.index = -1;
           end
        end
        
        function i = getSelectedIndex(obj)
            i = obj.current.index;
        end
        function unSetIndex(obj)
           obj.current.index = -1;
           obj.notify('indexChanged');
        end
        
        
        function it = getSelectedItemName(obj)
            list = obj.current.getList();
            it = list{obj.current.index};
        end
        
        function path = getPath(obj)
           path = cell(1,length(obj.stack) + 1);
           for i = 1:length(obj.stack)
              path{i} = obj.stack{i}.getLabel(); 
           end
           path{end} = obj.current.getLabel();
        end
        
        function selectItem(obj, index)
            item = obj.current.selectItem(index);
            
            if (~isempty(item))
                obj.stack{end+1} = obj.current;
                obj.current = item;
                obj.current.initIndex();
                obj.notify('listChanged');
            end
        end
        
        function goUpTo(obj, pathindex)
            len = length(obj.stack);
            if (pathindex <= len)
                for i = len:-1:pathindex
                    delete(obj.current);
                    obj.current = obj.stack{i};
                    obj.stack(i) = [];
                end
                obj.current.backOnTop();
            end
            obj.notify('listChanged');
        end
        
        function goUp(obj)
            if (~isempty(obj.stack) && obj.current.goUp())
                delete(obj.current);
                obj.current = obj.stack{end};
                obj.stack(end) = [];
                obj.notify('listChanged');
                obj.current.backOnTop();
            end
        end
        
        function c = getCurrent(obj)
            c = obj.current;
        end

        function str = getCurrentTooltipString(obj)
          str = obj.current.getToolTipString();
        end
        
        function keypress(obj, ev)
            if strcmp(ev.Key, 'return')
                if obj.current.index > 0
                    drawnow;
                    obj.selectItem(obj.current.index);
                end
            else
                obj.current.keypress(ev); 
            end
        end
        
        function shutDown(obj)
           obj.notify('shutdown'); 
        end
    end
    
    methods(Static)
        %{
        function im = constructInspectorModel(session , system, diagram , curve , displaysystems)
            
            im = InspectorModel(session);

            passon.preview = true;
            passon.session = session;            
            passon.killswitch = @() im.shutDown() ;
            passon.listreloaded = @() im.notify('listChanged');
            passon.im = im;
            passon.system = system;
            
            im.current = InspectorSystem(passon,  session.getSystemsPath() ,'Systems');
            im.current.index = -1;
            
            im.stack{1} = im.current;
            im.current = InspectorDiagram(passon, system.getPath() ,  system.getName());
            im.current.index = -1;
            
            if (~isempty(diagram))
                if(isobject(diagram)), diagram = diagram.diagrampath; end
                passon.diagram = choppath(diagram);
                im.stack{2} = im.current;
                im.current = InspectorCurve(passon, diagram ,  passon.diagram,[]);
                im.current.index = -1;
                
                if (~isempty(curve) && curve.hasFile() && ~curve.isEmpty())
                    passon.curvemanager = im.current.curvemanager;
                    passon.curve = curve;
                    im.stack{3} = im.current;
                    im.current = InspectorPoint(passon,curve , curve.getLabel());
                    im.current.index = -1;
                end
            end
            
        
        end
        %}
        
        %{
        function im = currentCurve(session , displaysystems) 
            if (nargin == 1), displaysystems = true; end
            im =   InspectorModel.constructInspectorModel(session, session.getSystem() , session.getCurveManager() , session.getCurrentCurve(), displaysystems);
            InspectorPanel.startWindow( session.getWindowManager() , im);
        end
        %}
        
        %{
        function im = currentDiagram(session, displaysystems)
            if (nargin == 1), displaysystems = true; end
            im = InspectorModel.constructInspectorModel(session, session.getSystem() , session.getCurveManager(), [] , displaysystems);
            InspectorPanel.startWindow( session.getWindowManager() , im);
        end
        
        %}
    end
end

%{
function str = choppath(path)
    if ((path(end) == '\') || (path(end) == '/')), path(end) = ''; end
    str = regexprep(path , '.*[/\\]' , '');
end
%}
