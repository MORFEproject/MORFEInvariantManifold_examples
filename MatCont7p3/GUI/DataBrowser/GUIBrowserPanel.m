classdef GUIBrowserPanel < handle
    
    properties
        panelhandle
        
        timestamp = 0;
        mainnode;
        previewpanel;
        resizetimer = [];
    end
    
    methods
        
        function obj = GUIBrowserPanel(parent, session , bm, showPath, varargin)
            set(parent,'Units' , 'Pixels');
            mainnode = LayoutNode(-1,-1,'vertical');
            obj.mainnode = mainnode;
            obj.panelhandle = uipanel(parent, 'Unit' , 'Pixels'  , varargin{:} ...
                , 'ResizeFcn' , @(o,e) obj.onResize(o,mainnode) , 'DeleteFcn' , @(o,e) obj.destructor());
            
            bpm = GUIBrowserPreviewModel(session, bm);
            
            if (showPath)
                mainnode.addGUIobject(1,1,   GUIBrowserPath( obj.panelhandle , bm ,  ... 
                    'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR) , 'halign' , 'l' , 'minsize' , [Inf,20]);
            end
            subnode = LayoutNode(12,1);
            subnode.addGUIobject(1,3, GUIBrowserList( obj.panelhandle , bm ,'BackgroundColor' , 'white') , 'minsize' , [Inf,Inf] , 'margin' , [5,2]);
            obj.previewpanel = GUIBrowserPreviewPanel(obj.panelhandle, bpm );
            subnode.addGUIobject(1,5, obj.previewpanel , 'minsize' , [Inf,Inf] , 'margin' , [5,2]);
            mainnode.addNode(subnode);
            mainnode.addGUIobject(1,1, GUIBrowserCommandPanel(obj.panelhandle, bm, 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR),'minsize' , [Inf,30]);

            bm.addlistener('shutdown' , @(o,e) close(parent) );
            mainnode.makeLayoutHappen( get( obj.panelhandle , 'Position') );
            obj.timestamp = tic();
            %set(obj.panelhandle, 'Units' , 'normalize');
            set(parent, 'ResizeFcn', @(o, e) obj.doResize());
            obj.resizetimer = timer('StartDelay', 0.2, 'TimerFcn', @(o, e) obj.syncPanel(get(parent, 'Position')));
            scrollAmmount = DefaultValues.LETTERDIMENSION(2) * 1;
            set(parent, 'Visible', 'on', 'WindowScrollWheelFcn', @(o, e) onScrollEvent(obj, e.VerticalScrollCount, scrollAmmount));
        end
        
        function doResize(obj)
            stop(obj.resizetimer); %reset timer.
            start(obj.resizetimer);
            
        end
        
        
        function syncPanel(obj, figureposition)
            set(obj.panelhandle, 'Position', [1 1 (figureposition(3:4)-1)]);
        end
        
        function onResize(obj,handle,mainnode)
            

            
            pos = get(handle, 'Position');
            if (~isempty(pos))
                mainnode.makeLayoutHappen( get(handle, 'Position'));
                obj.timestamp = tic();
            end
            
        end
        function destructor(obj)
           delete(obj.mainnode);
           delete(obj);
            
        end
        
        
    end
    methods(Static)
        function startWindow(wm , bm, varargin)

                f = wm.demandWindow('browser');

                set(f, 'Visible' , 'off'  );  %, 'WindowStyle' , 'modal');
                pos = get(f, 'Position');
                InspectorPanel(f , true , bm , 'Position' , [0 0 pos(3) pos(4)] , varargin{:});
                set(f, 'Visible' , 'on');
        end
        function startWindowWoPath(wm  , bm, windowlabel ,  varargin)

                f = wm.demandWindow(windowlabel);

                set(f, 'Visible' , 'off'  );  %, 'WindowStyle' , 'modal');
                pos = get(f, 'Position');
                InspectorPanel(f , false , bm , 'Position' , [0 0 pos(3) pos(4)] , varargin{:});
                set(f, 'Visible' , 'on');
        end        
    end
    
end
function onScrollEvent(panel, direction, amount)
    %direction: 1 if scroll down, -1 if scroll up.
    if isempty(panel.previewpanel.slider) || ~isvalid(panel.previewpanel.slider); return; end
    
    slider = panel.previewpanel.slider.handle;
    newvalue = slider.Value - direction*amount;
    newvalue = min(max(slider.Min, newvalue), slider.Max);
    slider.Value = newvalue;
    slider.Callback(slider, []);
    

   
end

