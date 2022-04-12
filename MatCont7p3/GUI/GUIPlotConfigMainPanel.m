classdef GUIPlotConfigMainPanel
    
    properties
        panelhandle;
        
        subpanel
        mainnode = [];
    end
    
    
    methods
        function obj = GUIPlotConfigMainPanel(wm, plotconf, plotposition ,varargin)
            windowlabel = ['plotlayout' num2str(plotconf.getDimension()) 'd'];
            
            parent = wm.demandWindow(windowlabel);
            
            
            
            obj.panelhandle = uipanel(parent, 'Unit' , 'Pixel'   , varargin{:});
            
            mainnode = LayoutNode(-1,-1,'vertical');
            
            obj.mainnode = mainnode;
            
            for i = 1:plotconf.getDimension()
                obj.subpanel{i} =  GUIPlotConfigPanel(obj.panelhandle , plotconf , i);
                mainnode.addHandle(4,0,obj.subpanel{i}.panelhandle , 'minsize' , [Inf Inf]);
            end
            
            
            
            subnode = LayoutNode(1,0);
            
            subnode.addHandle(0,1,uicontrol(obj.panelhandle,'Style','pushbutton' , 'String' , 'OK',...
                'Callback' , @(o,e) closeOK(parent, plotconf) ), 'minsize' , [100 , 25],'margin' , [10 2]);
            mainnode.addNode(subnode);
            
            set(obj.panelhandle,  'ResizeFcn' , @(o,e) obj.doLayout());
            set(obj.panelhandle, 'Units' , 'normalize');
            obj.doLayout();
            
            set(obj.panelhandle, 'DeleteFcn' , @(o,e) obj.destructor());
            
            
            position = get(parent, 'Position');
            newposition = [plotposition(1:2), position(3:4)];
            newposition(2) = newposition(2) + round(max(plotposition(4)/2 - position(4)/2, 0));
            newposition(1) = newposition(1) + round(max(plotposition(3)/2 - position(3)/2, 0));
            set(parent, 'Visible' , 'on', 'Position', newposition);
        end
        
        function waitForClose(obj)
            waitfor(obj.panelhandle);
        end
        
        function doLayout(obj)
            units = get(obj.panelhandle, 'Units');
            
            set(obj.panelhandle,'Units', 'Pixels');
            obj.mainnode.makeLayoutHappen( get(obj.panelhandle, 'Position'));
            set(obj.panelhandle,'Units',units);
        end
        
        function destructor(obj)
            delete(obj.mainnode);
        end
        
    end
end

function closeOK(fighandle, plotconf)
    delete(fighandle);
    %disp(plotconf);
    %cla(plotconf.axeshandle);  clear plot
end