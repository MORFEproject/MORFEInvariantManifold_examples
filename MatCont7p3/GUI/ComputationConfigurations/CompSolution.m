classdef CompSolution
    properties(Hidden)
       settings; 
       compbranch;
       name;
    end
    
    
    methods
        function obj = CompSolution(settings, compbranch) 
            obj.settings = settings.copy();
            obj.compbranch = compbranch;
            obj.name = '<no name set>';
        end
        
        function list = listSwitches(obj, verbose)
           list = [];
        end

        
        function plotdim(obj, dim)
            oi = obj.compbranch.getOutputInterpreter();
            [~, ~, plotselections] = oi.interpret(obj, obj.settings, obj.compbranch);
            
            
            wm = GUIWindowManager([]);
            
            previousplots = obj.settings.previousplots;
            plotconf = previousplots.recover(dim, wm);
            if isempty(plotconf)
                plotconf = GUIPlotConf(dim, wm);
                
            end
            
            plotconf.configurePlotSelection(plotselections);
            plotconf.setSolutionHandler(CLSolutionHandler.fromSingleSolution(obj));
            plotconf.generateFigure();
            
            layoutpanel = plotconf.showLayoutWindow();
            
            layoutpanel.waitForClose();
            plotconf.isValid();
            out = plotconf.getOutputter(oi.getPlotter(), []);
            if ~isempty(out)
                out.plotSolution(obj);
            end
        end
        

        
        function plot(obj)
            obj.plotdim(2);
        end
        
        function plot3(obj)
           obj.plotdim(3); 
        end
        
        function n = getName(obj)
            n = obj.name;
        end
        
        function save(obj, filename) 
        end
        
        function [b, msg] = showTable(obj, windowname, session)
            b = 0; msg = 'not implemented';
        end
        
        function restoreGlobals(obj)
        end
        
        
    end
    
    methods(Static)
         function obj = load(filename)
           data = load(filename); 
           
           if ~isfield(data, 'gui')
                obj = LegacySolutionLoader(filename);
           else
                obj = data.gui.loader(data);
           end
        end       
    end
end