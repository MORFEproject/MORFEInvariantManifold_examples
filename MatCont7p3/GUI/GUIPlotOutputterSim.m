classdef GUIPlotOutputterSim < handle
    properties
        dimension;
        
        fx
        fy
        fz
        
        axeshandle
        group 
        getData
        
        plotconfig;
        pointloader;
        curvelabel;       
        
        previouspoint;
    end
    
    methods
        function obj = GUIPlotOutputterSim(axeshandle, fx, fy, fz, ~, plotconfig, pointloader, curvelabel)
            obj.fx = fx;
            obj.fy = fy;
            obj.fz = fz;
            obj.axeshandle = axeshandle;
            
            if isempty(fz)
                obj.dimension = 2;
                obj.getData = @(data, s, ind) {fx(data{:}, s, ind, ind(end)),  fy(data{:}, s, ind, ind(end))};

            else
                obj.dimension = 3;
                obj.getData = @(data, s, ind) {fx(data{:}, s, ind, ind(end)),  fy(data{:}, s, ind, ind(end)), fz(data{:}, s, ind, ind(end))};
            end  
            

    
            obj.pointloader = pointloader;
            obj.plotconfig = plotconfig; 
            obj.curvelabel = curvelabel;
            
            obj.group = hggroup('Parent' , obj.axeshandle );
            obj.previouspoint = {[], [], []};
            
        end
        
        function data = addPreviousPoint(obj, data)
            for i = 1:obj.dimension
                packet = data{i};
                data{i} = [obj.previouspoint{i}; packet(:)];
            end
            
        end
        function output(obj, data, s, ind)

            data = obj.getData(data, s, ind);
            data = obj.addPreviousPoint(data);
            obj.previouspoint = cellfun(@(x) {x(end)}, data);
            
            line(obj.axeshandle, data{:}, 'Parent', obj.group,  obj.plotconfig.curve{:});
            
            
        end
 
    
        function plotSolution(obj, solution)
           obj.output({solution.t, solution.y}, [], 1:length(solution.t));
           
           for k = 1:length(solution.tE)
                obj.outputPoint(solution.tE(k), solution.yE(k,:), ['E' num2str(solution.iE(k))]);
           end
        end
        
        function outputPoint(obj, t, y, label)
            labeldata = obj.getData({t, y}, [], 1);
            line(labeldata{:},'Parent', obj.group, obj.plotconfig.specialpoint{:});
            text(labeldata{:}, label, 'Parent' ,obj.axeshandle,  'ButtonDownFcn', ...
                    @(o,e) selectPoint(o, obj.pointloader, obj.curvelabel), 'UserData', {t, y, label})
        
        end
        
        function vector = getSkewDirection(obj)
            d = axis(obj.axeshandle);
            
            vector = [0 0 0];
            for dim = 1:obj.dimension
                i = dim - 1;
                vector(dim) = d(2 + 2*i) - d(1 + 2*i);
            end
            vector = 0.012*vector;            
        end
        
        
    end
end

function selectPoint(guiobject, pointloader, label)
    if isempty(pointloader)
        disp(repmat('-', 1, 80))
        disp(guiobject.UserData);
        if isfield(guiobject.UserData, 'data')
           disp(guiobject.UserData.data); 
        end
        disp(repmat('-', 1, 80))
    else
        if isvalid(pointloader)
            pointloader.loadPointSim(label, guiobject.UserData);
        else
            fprintf(2, 'Connection with current session has been lost, unable to load in point\n');
            
        end
    end

end



