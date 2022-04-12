%default output interpreter for continuation curves
classdef CLContOutputInterpreter
    
    methods
        function [coord, param, period, meta] = getPoint(~, solution, index)
            meta = struct();
            settings = solution.settings;
            parameters = settings.getSetting('parameters');
            param = parameters.getValue();
            ap = parameters.getActive();
            dim = length(settings.coord);
            coord = solution.x(1:dim, index);
            period = [];
            
            hasperiod = solution.settings.getVal('Period', 0) || solution.compbranch.hasPeriod();
            
            [ntst, ncol] = solution.getDiscretization();
            leap = hasperiod; %leap over state space data
            leap = leap + (ntst*ncol + 1)*dim;
            
            if hasperiod
                period = solution.x(leap, index);
            end
            
            param(ap) = solution.x(leap+1:leap+length(ap), index);
        end
        
        function plotter = getPlotter(~)
            plotter = @GUIPlotOutputter;
        end
        function renderer = getPreviewPanelRenderer(~)
            renderer = @(sol, pp, pm) GUIPreviewRendererSolutions.renderContinuationCurve(sol, pp, pm);
        end
        
        
        
        
        function [omap, numconfig, plotsel] = interpret(obj, solution, settings, compbranch)
            
            numconfig = settings.numericconfig;
            numconfig.reset();
            plotsel = PlotSelections();
            
            if isempty(solution)
                omap = ContOutputMap(0, 0, 0, 0);
            else
                omap = ContOutputMap(size(solution.x, 1), size(solution.v, 1), size(solution.h, 1), size(solution.f, 1));
            end
            
            system = settings.system;
            coordinates = system.getCoordinates();
            dim = length(coordinates);
            parameters = settings.getSetting('parameters');
            
            if isempty(solution)
                %fprintf(2, 'WARNING: ntst/ncol should be extracted from lds, not settings\n\n'); FIXME TODO
                ntst = settings.getVal('ntst', 0);
                ncol = settings.getVal('ncol', 0);
            else
                [ntst, ncol] = solution.getDiscretization();
            end
            
            
            omap = obj.doX(settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel);
            omap = obj.doH(settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel);
            omap = obj.doF(settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel);
            
            
            numconfig.declareCategory('npoints', 'NPoints',1e7, false)
            plotsel.declareCategory('npoints', 'NPoints',1e7, true)
            numconfig.setLabels('npoints', {'Npoints', @(x, h, f, s, ind, i) i});
            plotsel.declareSubCategory('npoints', 'Npoints', @(x, h, f, s, ind, i) ind);
            
            plotsel.setSolutionLabel(compbranch);
            numconfig.configDone();
        end
        
        function [index, omap] = doCycle(obj, coordinates, ntst, ncol, dim, omap, numconfig, plotsel)
            
            for i = 1:length(coordinates)
                omap.x{i} =  coordinates{i};  %['coordinate ' coordinates{i}];
            end
            
            
            
            if ntst == 0
                numconfig.declareCategory('coordinates', 'Coordinates',1e0, true);
                plotsel.declareCategory('coordinates', 'Coordinates', 1e0);
                
                labels = cell(length(coordinates), 2);
                for index = 1:length(coordinates)
                    labels(index, :) = {coordinates{index}, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('coordinates', coordinates{index}, @(x, h, f, s, ind, i) x(index, ind))
                end
                numconfig.setLabels('coordinates', labels);
                
            else
                plotsel.declareCategory('coordinates', 'Coordinates', 1e0);
                for index = 1:length(coordinates)
                    coordinate = coordinates{index};
                    plotsel.declareSubCategory('coordinates', coordinate);
                    
                    selection = (0:ntst*ncol)*dim + index;
                    plotsel.declareItem('coordinates', coordinate, 'Default', coordinate, @(x, h, f, s, ind, i) x(selection, ind), true);
                    plotsel.declareItem('coordinates', coordinate, 'Min', sprintf('min(%s)', coordinate), @(x, h, f, s, ind, i) min(x(selection, ind)));
                    plotsel.declareItem('coordinates', coordinate, 'Max', sprintf('max(%s)', coordinate), @(x, h, f, s, ind, i) max(x(selection, ind)));
                end
                
                
                
            end
            
            
            omap.x(dim+1:dim+(ntst*ncol)*dim) = repmat(coordinates, 1, ntst*ncol);
            index = dim*(ntst*ncol+1) + 1;
            
            
        end
        
        
        function [index, omap] = doParameters(obj, index, settings, compbranch, parameters, omap, numconfig, plotsel)
            
            if isempty(parameters); return;  end
            
            actives = parameters.getActive();
            numconfig.declareCategory('parameters', 'Parameters',1e1, true)
            plotsel.declareCategory('parameters', 'Parameters',1e1)
            labels = cell(length(actives), 2);
            parameternames = parameters.parameters;
            parametervalues = parameters.getValue();
            
            j = 1;
            for ii = 1:length(parameternames)
                parametername = parameters.parameters{ii};
                
                if all(ii ~= actives)
                    plotsel.declareSubCategory('parameters', parametername, @(x, h, f, s, ind, i) repmat(parametervalues(ii), 1, length(ind)));
                else
                    omap.x{index} = parametername; %['param ' parameters.parameters{active}];
                    labels(j, :) = {parametername, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('parameters', parametername, @(x, h, f, s, ind, i) x(index, ind))
                    index = index + 1;
                    j = j + 1;
                end
                
            end
            numconfig.setLabels('parameters', labels);
            
        end
        
        function omap = doX(obj, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            [index, omap] = doCycle(obj, coordinates, ntst, ncol, dim, omap, numconfig, plotsel);
            
            period = settings.getVal('Period', 0) || compbranch.hasPeriod();
            if period
                omap.x{index} = 'Period';
                numconfig.declareCategory('period', 'Period',1e2, true)
                plotsel.declareCategory('period', 'Period',1e2, true)
                
                plotsel.declareSubCategory('period', 'Period', @(x, h, f, s, ind, i) x(index, ind))
                numconfig.setLabels('period', {'Period', @(x, h, f, s, ind, i) x(index, i)});
                index = index + 1;
            end
            
            [~, omap] = obj.doParameters(index, settings, compbranch, parameters, omap, numconfig, plotsel);
        end
        
        function omap = doH(obj, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            omap.h{1} = 'stepsize';
            numconfig.declareCategory('stepsize', 'Stepsize',1e5, false)
            plotsel.declareCategory('stepsize', 'Stepsize',1e5, true)
            numconfig.setLabels('stepsize', {'Stepsize', @(x, h, f, s, ind, i) h(1, i)});
            plotsel.declareSubCategory('stepsize', 'stepsize', @(x, h, f, s, ind, i) h(1, ind))
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            omap.h{2} = 'correction';
            %numconfig.declareCategory('corrections', 'Corrections',1e5+1, false)
            %numconfig.setLabels('corrections', {'Corrections', @(x, h, f, s, ind, i) h(2, i)});
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            uf = settings.getVal('userfunctions', []);
            index = 2;
            if ~isempty(uf) && ~isempty(find([uf.ufdata.state]))
                numconfig.declareCategory('uf', 'User Functions',1e6, false)
                plotsel.declareCategory('uf', 'User Functions',1e6)
                ufdata = uf.ufdata;
                actives = find([ufdata.state]);
                labels = cell(length(actives), 2);
                j = 1;
                
                for k = 1:length(ufdata)
                    index = index + 1;
                    if any(k == actives)
                        omap.h{index} = strip(ufdata(k).label);
                        labels(j, :) = {strip(ufdata(k).label), @(x, h, f, s, ind, i) h(index, i)};
                        j = j + 1;
                        plotsel.declareSubCategory('uf',strip(ufdata(k).label) , @(x, h, f, s, ind, i) h(index, ind))
                    else
                        omap.h{index} = '--';
                    end
                end
                
                numconfig.setLabels('uf', labels);
            end
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [index, omap] = obj.doTests(index, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel);
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %branch parameters
            if isempty(parameters); return;  end
            bp = parameters.getBranch();
            if ~isempty(bp)
                numlabels = cell(0, 2);
                numconfig.declareCategory('bptf', 'Branch Parameter Test Functions',1001, false)
                plotsel.declareCategory('bptf', 'Branch Parameter Test Functions', 1001)
                for bpi = bp
                    index = index + 1;
                    plotsel.declareSubCategory('bptf', parameters.parameters{bpi}, @(x, h, f, s, ind, i) h(index, ind));
                    numlabels(end+1, :) = {parameters.parameters{bpi}, @(x, h, f, s, ind, i) h(index, i)};
                    omap.h{index} = ['bp ' parameters.parameters{bpi}];
                end
                numconfig.setLabels('bptf', numlabels);
                
            end
        end
        
        function [index, omap] = doTests(obj,index, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            [labels, internallabels] = compbranch.getTestLabels();
            
            numlabels = cell(0, 2);
            plotlabels = cell(0, 2);
            for i = 1:length(labels)
                if settings.getValue(internallabels{i})
                    if strcmp(labels{i}, 'BPC')  %LC case
                        for j = 1:4
                            index = index + 1;
                            omap.h{index} = [labels{i} num2str(j)];
                            numlabels(end+1, :) = {[labels{i} num2str(j)], @(x, h, f, s, ind, i) h(index, i)};
                            plotlabels(end+1, :) = {[labels{i} num2str(j)], @(x, h, f, s, ind, i) h(index, ind)};
                        end
                    else
                        index = index + 1;
                        omap.h{index} = labels{i};
                        numlabels(end+1, :) = {labels{i}, @(x, h, f, s, ind, i) h(index, i)};
                        plotlabels(end+1, :)= {labels{i}, @(x, h, f, s, ind, i) h(index, ind)};
                    end
                end
            end
            if ~isempty(numlabels)
                numconfig.declareCategory('tf', 'Test Functions',1e3, false)
                plotsel.declareCategory('tf', 'Test Functions',1e3)
                numconfig.setLabels('tf', numlabels);
                for k = 1:size(plotlabels, 1)
                    plotsel.declareSubCategory('tf', plotlabels{k, 1}, plotlabels{k, 2});
                end
            end
        end
        
        function omap = doF(obj, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            
            index = 1;
            if ntst > 0
                omap.f(1:ntst+1) = repmat({'(mesh)'}, 1, ntst+1);
                index = ntst+1+1;
                
                PRC = settings.getVal('PRCenabled', 0);
                dPRC = settings.getVal('dPRCenabled', 0);
                
                if PRC
                    omap.f(index:index + (ntst*ncol + 0)) = repmat({'PRC'}, 1, ncol*ntst+1);
                    index = index + (ntst*ncol + 1);
                end
                if dPRC
                    omap.f(index:index + (ntst*ncol + 0)) = repmat({'dPRC'}, 1, ncol*ntst+1);
                    index = index + (ntst*ncol + 1);
                end
                
            end
            
            if settings.getVal('eigenvalues', 0)
                omap.f(index:index+dim-1) = repmat({'eig'}, 1, dim);
                numconfig.declareCategory('eig', 'Eigenvalues',1e4, false)
                plotsel.declareCategory('eig', 'Eigenvalues',1e4)
                
                labels = cell(dim*2, 2);
                for k = 1:dim
                    plabel = sprintf('eig%i', k);
                    plotsel.declareSubCategory('eig', plabel);
                    labels(k, :) =  {sprintf('Re[%i]', k), @(x, h, f, s, ind, i) real(f(index, i))};
                    plotsel.declareItem('eig', plabel, 'Re', sprintf('Re[%i]', k), @(x, h, f, s, ind, i) real(f(index, ind)));
                    labels(k + dim, :)   =  {sprintf('Im[%i]', k), @(x, h, f, s, ind, i) imag(f(index, i))};
                    plotsel.declareItem('eig', plabel, 'Im', sprintf('Im[%i]', k), @(x, h, f, s, ind, i) imag(f(index, ind)));
                    index = index + 1;
                end
                
                numconfig.setLabels('eig', labels);
                
                
            elseif settings.getVal('multipliers', 0)
                numconfig.declareCategory('mult', 'Multipliers',1e4+1, false)
                plotsel.declareCategory('mult', 'Multipliers',1e4+1)
                omap.f(index:index+dim-1) = repmat({'mult'}, 1, dim);
                
                labels = cell(dim*2, 2);
                for k = 1:dim
                    plabel = sprintf('mult%i', k);
                    plotsel.declareSubCategory('mult', plabel);
                    labels(k, :) =  {sprintf('Mod[%i]', k), @(x, h, f, s, ind, i) abs(f(index, i))};
                    plotsel.declareItem('mult', plabel, 'Mod', sprintf('Mod[%i]', k), @(x, h, f, s, ind, i) abs(f(index, ind)));
                    labels(k + dim, :)   =  {sprintf('Arg[%i]', k), @(x, h, f, s, ind, i) rad2deg(angle(f(index, i)))};
                    plotsel.declareItem('mult', plabel, 'Arg', sprintf('Arg[%i]', k), @(x, h, f, s, ind, i) rad2deg(angle(f(index, ind))));
                    index = index + 1;
                end
                
                numconfig.setLabels('mult', labels);
                
            end
        end
        
    end
    
end
