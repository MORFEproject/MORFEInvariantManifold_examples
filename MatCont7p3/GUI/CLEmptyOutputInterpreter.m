%Empty output interpreter, for testing only.

classdef CLEmptyOutputInterpreter
    
    methods
        function [coord, param, period, meta] = getPoint(~, solution, index)
            assert(false, 'refuse to execute -- Testing only');
        end
        
        function plotter = getPlotter(~)
            assert(false, 'refuse to execute -- Testing only');
        end
        function renderer = getPreviewPanelRenderer(~)
            assert(false, 'refuse to execute -- Testing only');
        end
        
        
        function [omap, numconfig, plotsel] = interpret(obj, solution, settings, compbranch)
            fprintf(2, 'empty interpretation call \n');           
            numconfig = settings.numericconfig;
            numconfig.reset();
            plotsel = PlotSelections();
            
            if isempty(solution)
                omap = ContOutputMap(0, 0, 0, 0);
            else
                omap = ContOutputMap(size(solution.x, 1), size(solution.v, 1), size(solution.h, 1), size(solution.f, 1));
            end
        
        end
    end
end
