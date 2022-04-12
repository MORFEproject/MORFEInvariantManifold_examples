classdef GUISimCurveTable

    methods(Static)
        function previewPanel(c, name, session)
            f = figure('Visible' , 'off' , 'NumberTitle' , 'off' , 'Name' , name, 'tag', 'matcont');
            if ~isempty(session)
                listener = session.addlistener('shutdown', @(o, e) close(f));
            else
                listener = [];
            end
            mainnode = LayoutNode(-1,-1,'vertical');
            
            mainnode.setOptions('Add' , true);
            h = GUISimCurveTable.makeTable(f, c.t, c.y, c.settings.system);
            
            mainnode.addHandle(10 , 1 , h, 'minsize' , [Inf, Inf]);
            mainnode.makeLayoutHappen( get(f,'Position'));
            set(f,'Visible' , 'on' , 'ResizeFcn' , @(o,e)  GUISimCurveTable.layoutResize( o, e), 'DeleteFcn' , @(o,e) deleteobjects(mainnode, listener), 'UserData', timer('StartDelay', 0.2, 'TimerFcn', @(~, ~) mainnode.makeLayoutHappen( get(f,'Position')) ));
        end
        
        function layoutResize(o, e)
            stop(o.UserData);
            start(o.UserData);
            
        end
        
        function handle = makeTable(parent, t, y, system, varargin)
            rownames = [system.getTimeName(), system.getCoordinates()];
            maxlen = max(cellfun(@(x) length(x), rownames));
            rownames = cellfun(@(x) {pad(x, maxlen)}, rownames);
            
           
            handle = uitable(parent, 'Data' , num2cell([t';y']) , 'RowName', rownames , varargin{:});
        end
        
    end
end


function deleteobjects(o1, o2)
    delete(o1);
    delete(o2);

end
