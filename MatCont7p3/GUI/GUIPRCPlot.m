classdef GUIPRCPlot
   
    properties
       ax 
       names;
       color;
       lightcolor;
       derivative;
       
       solutionhandle;
    end
    
    methods
        function obj = GUIPRCPlot(session, fig, names, color, lightcolor, derivative)
           if nargin == 0
              fig = figure();  
           end
           
            
            set(fig, 'MenuBar', 'figure', 'color',[1 1 1], 'UserData', [], 'KeyPressFcn', @(o,e) SessionOutputKeypress(o, e));
            
            obj.ax = axes('Parent' , fig);
                      
            
            obj.names = names;
            obj.color = color;
            obj.lightcolor = lightcolor;
            obj.derivative = derivative;
            
            if ~isempty(session)
                obj.solutionhandle = session.solutionhandle;
                
                mhandle = uimenu(fig, 'Label' , 'MatCont');
                uimenu(mhandle , 'Label' , 'Redraw Current Curve' , 'Callback' , @(o,e) obj.plotCurrentCurve());
                
            end
            
            
            
            obj.reset();   
        end
        
        function b = doInit(obj)
           obj.reset();
           b = 1; 
        end
        function reset(obj)
            d = axis(obj.ax);
            p = get(obj.ax, 'View');
            xmode = get(obj.ax,'Xdir');
            ymode = get(obj.ax,'YDir');
            cla(obj.ax)
            obj.ax.UserData = [];
            set(obj.ax,'View',p,'XDir',xmode,'YDir',ymode);
            axis(obj.ax, d);
             title(obj.ax, obj.names{1});
            xlabel(obj.ax, obj.names{2});
            ylabel(obj.ax, obj.names{3});   
            line(obj.ax, [0 1], [0 0], 'Color' ,'k','LineWidth',1);
            axis(obj.ax, 'tight');
             
        end
        
        

        

        
        function plotPRC(obj, finemesh, data)
            obj.ax.UserData = line(obj.ax, finemesh, data, 'Color',  obj.color, 'LineWidth', 1);
        end
        
        function output(obj, data, s, ind)
            
            global lds;
            start = lds.ntst + 1;
            nrfinemesh = lds.ntst*lds.ncol + 1;
            
            if lds.CalcPRC && ~obj.derivative
                for i = ind
                    set(obj.ax.UserData, 'color', obj.lightcolor);
                    obj.plotPRC(lds.finemsh, data{3}(start+1:start+nrfinemesh, i))
                end
                
            end
            start = start + nrfinemesh;
            if lds.CalcdPRC && obj.derivative
                for i = ind
                    set(obj.ax.UserData, 'color', obj.lightcolor);
                    obj.plotPRC(lds.finemsh, data{3}(start+1:start+nrfinemesh, i))
                end
            end
        end
        
        function plotCurrentCurve(obj)
            obj.reset();
            if isempty(obj.solutionhandle); return; end
            if isempty(obj.solutionhandle.solution); return; end
            try
                solution = obj.solutionhandle.solution;
                lds = solution.globals.lds; %assume continuation curve.
                start = lds.ntst + 1;
                nrfinemesh = lds.ntst*lds.ncol + 1;
                if lds.CalcPRC && ~obj.derivative
                    for i = 1:size(solution.f, 2)
                        obj.plotPRC(lds.finemsh, solution.f(start+1:start+nrfinemesh, i))
                    end
                    
                end
                start = start + nrfinemesh;
                if lds.CalcdPRC && obj.derivative
                    for i = 1:size(solution.f, 2)
                        obj.plotPRC(lds.finemsh, solution.f(start+1:start+nrfinemesh, i))
                    end
                end
            catch
            end
        
        end
    end
    
end



