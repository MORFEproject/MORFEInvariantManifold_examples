%default output interpreter for continuation curves
classdef CLSimOutputInterpreter
    
    methods
        
        %{
        %unused.
        function [coord, param, period, meta] = getPoint(~, solution, index)
        end
        %}
        
        function plotter = getPlotter(~)
           plotter = @GUIPlotOutputterSim;
        end
        
        
        function renderer = getPreviewPanelRenderer(~)
           renderer = @(sol, pp, pm) GUIPreviewRendererSolutions.renderSimulation(sol, pp, pm);
        end
        
        function [omap, numconfig, plotsel] = interpret(~, solution, settings, compbranch)
            omap = []; %used for continuation curves.

            numconfig = settings.numericconfig;
            numconfig.reset();
            plotsel = PlotSelections();
            
            system = settings.system;
            coordinates = system.getCoordinates();
            dim = length(coordinates);
            parameters = system.getParameters();
            parametervalues = settings.getValue('parameters');
            
            numconfig.declareCategory('coordinates', 'Coordinates',1e0, true);
            plotsel.declareCategory('coordinates', 'Coordinates', 1e0);
            
            labels = cell(dim, 2);
            for index = 1:dim
                labels(index, :) = {coordinates{index}, @(t, y, s, ind, i) y(i, index)};
                plotsel.declareSubCategory('coordinates', coordinates{index}, @(t, y, s, ind, i) y(ind, index))
            end
            numconfig.setLabels('coordinates', labels);
            
            
            numconfig.declareCategory('parameters', 'Parameters',1e2, false)
            plotsel.declareCategory('parameters', 'Parameters',1e2)
            labels = cell(length(parameters), 2);
            for index = 1:length(parameters)
                labels(index, :) = {parameters{index}, @(t, y, s, ind, i) parametervalues(index) };
                plotsel.declareSubCategory('parameters', parameters{index},  @(t, y, s, ind, i) repmat(parametervalues(index), 1, length(ind)));
            end
            numconfig.setLabels('parameters', labels);
            
            timevar = system.getTimeName();
            numconfig.declareCategory('time', 'Time',101, true);
            numconfig.setLabels('time', {timevar,  @(t, y, s, ind, i) t(i)});
            plotsel.declareCategory('time', 'Time', 101);
            plotsel.declareSubCategory('time', timevar, @(t, y, s, ind, i) t(ind));
            
           
            if compbranch.hasTarget()
                targetsetting = settings.getSetting('coord_target');
                if ~isempty(targetsetting) && targetsetting.isVisible()
                    target = targetsetting.getValue();
                else
                    target = settings.coord;
                end
                target = target(:)';
                name = compbranch.hasTarget();
                numconfig.declareCategory(name, name,1e4, true);
                numconfig.setLabels(name, {name,  @(t, y, s, ind, i) norm(y(i,:) - target)});
            end
            
            ufdata = system.ufdata;
            if ~isempty(ufdata) && ~isempty(find([ufdata.valid]))
                paramcell = num2cell(parametervalues);
                handles = system.handle(); %evaluate handles
                %out{1} = @init;out{2} = @fun_eval;out{3} = @jacobian;out{4} = @jacobiap;out{5} = @hessians;%
                %out{6} = @hessiansp;out{7} = @der3;out{8} = [];out{9} = [];
                %out{10}= @Userfun1;out{11}= @Userfun2; ... ...             
                handles(1:9) = [];
                numconfig.declareCategory('uf', 'User Functions',1e6, false)
                plotsel.declareCategory('uf', 'User Functions',1e6)
                actives = find([ufdata.valid]);
                labels = cell(length(actives), 2);
                j = 1;
                
                for k = 1:length(ufdata)
                    index = index + 1;
                    if any(k == actives)
                        userfun = handles{k};
                        labels(j, :) = {strip(ufdata(k).label), @(t, y, s, ind, i) userfun(t(i), y(i, :)', paramcell{:})};
                        j = j + 1;
                        plotsel.declareSubCategory('uf',strip(ufdata(k).label) , @(t, y, s, ind, i)  arrayfun(@(z)userfun(t(z), y(z, :)', paramcell{:}), ind));
                    end
                end
                
                numconfig.setLabels('uf', labels);
            end
            
            
            
            plotsel.setSolutionLabel(compbranch);
            numconfig.configDone();
        end
        
    end
end
