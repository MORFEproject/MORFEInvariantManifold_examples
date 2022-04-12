%output interpreter for connecting orbit curves
classdef CLContOutputInterpreterConnection < CLContOutputInterpreter
    
    
    properties
       heteroclinic; 
        
    end
    
    methods
        
        function obj = CLContOutputInterpreterConnection(varargin)
           obj.heteroclinic = nargin ~= 0; 
        end
        
        function [coord, param, period, meta] = getPoint(obj, solution, index)
            
            settings = solution.settings;
            parameters = settings.getSetting('parameters');
            param = parameters.getValue();
            ap = parameters.getActive();
            dim = length(settings.coord);
            [ntst, ncol] = solution.getDiscretization();

            jump = (ncol*ntst + 1)*dim;
            coord = solution.x(jump+1:jump+dim,index);
            jump = jump + dim + 1;
            
            
            if obj.heteroclinic
                coord_target = solution.x(jump:jump+dim-1,index);
                jump = jump + dim;
            else
                coord_target = [];
            end
            param(ap) = solution.x(jump:jump+length(ap)-1, index);
            
            
            cparam = settings.getSetting('homoclinicparameters');
            homactive = cparam.getActive();
            %check if T (period) is selected in extravec/homoclinic
            %parameters, T is on index '1', eps0 on 2, eps1 on 3
            if any(homactive == 1)
                period = solution.x(jump+length(ap), index);
                
            else
                %if period not active, take period from initial point.
                period = solution.s(1).data.T;
                
            end
            
            eps0 = cparam.value(2);
            eps1 = cparam.value(3);

            i0 = find(homactive == 2);
            if ~isempty(i0)
                eps0 = solution.x(jump+length(ap) - 1 + i0, index);
            end
            i1 = find(homactive == 3);
            if ~isempty(i1)
                eps1 = solution.x(jump+length(ap) - 1 + i1, index);
            end
            
            meta = struct('eps0', eps0, 'eps1', eps1, 'coord_target', coord_target);
        end
        

        %{
        function renderer = getPreviewPanelRenderer(~)
            renderer = @(sol, pp, pm) GUIPreviewRendererSolutions.renderContinuationCurve(sol, pp, pm);
        end
        %}
        

        function omap = doX(obj, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            [index, omap] = doCycle(obj, coordinates, ntst, ncol, dim, omap, numconfig, plotsel);
            
            numconfig.declareCategory('coordinates_saddle_start', 'Coordinates (saddle)',1e1+1, false);
            plotsel.declareCategory('coordinates_saddle_start', 'Coordinates (saddle)', 1e1+1);
            labels = cell(length(coordinates), 2);
            labelindex = 1;
            for j = 1:length(coordinates)
                omap.x{index} = ['start ' coordinates{j}];
                labels(labelindex, :) = {coordinates{j}, @(x, h, f, s, ind, i) x(index, i)};
                plotsel.declareSubCategory('coordinates_saddle_start', coordinates{j}, @(x, h, f, s, ind, i) x(index, ind));
                index = index + 1;
                labelindex = labelindex + 1;
            end
            numconfig.setLabels('coordinates_saddle_start', labels);
            
            
            if obj.heteroclinic
                numconfig.declareCategory('coordinates_saddle_end', 'Coordinates (target saddle)',1e1+1, false);
                plotsel.declareCategory('coordinates_saddle_end', 'Coordinates (target saddle)', 1e1+1);
                labels = cell(length(coordinates), 2);
                labelindex = 1;
                for j = 1:length(coordinates)
                    omap.x{index} = ['end ' coordinates{j}];
                    labels(labelindex, :) = {coordinates{j}, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('coordinates_saddle_end', coordinates{j}, @(x, h, f, s, ind, i) x(index, ind));
                    index = index + 1;
                    labelindex = labelindex + 1;
                end
                numconfig.setLabels('coordinates_saddle_end', labels);
                
            end
            [index, omap] = obj.doParameters(index, settings, compbranch, parameters, omap, numconfig, plotsel);
            
            
                    
            cparam = settings.getSetting('homoclinicparameters');
            if isempty(cparam); return; end
            active = cparam.getActive();
            if active
                numconfig.declareCategory('conorb_parameters', 'Connection Parameters',1e1+2, true);
                plotsel.declareCategory('conorb_parameters', 'Connection Parameters', 1e1+2);
                labels = cell(length(active), 2);
                labelindex = 1;
                for j = active
                    paramname = cparam.parameters{j};
                    omap.x{index} = paramname;
                    
                    labels(labelindex, :) = {paramname, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('conorb_parameters', paramname, @(x, h, f, s, ind, i) x(index, ind));
                    labelindex = labelindex + 1;
                    index = index + 1;
                end
                numconfig.setLabels('conorb_parameters', labels);
                
            end
            
        end
        
        
        function [index, omap] = doTests(obj,index, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            
            [labels, internallabels] = compbranch.getTestLabels();
            
            numlabels = cell(0, 2);
            plotlabels = cell(0, 2);
            
            fs = homoclinic; singmat = feval(fs{9});
            mask = cellfun(@(x) settings.getValue(x)~=0,  internallabels, 'UniformOutput', 0);
            mask = [mask{:}];
            labels = labels(mask);
            internallabels = internallabels(mask);
            
            singmat = singmat(mask, :);
            mask = sum(singmat == 0) ~= 0;
            singmat = singmat(:, mask);
            
            
            baseindex = index;
            
            for i = 1:length(labels)
                label = labels{i};
                
                testfi = find(singmat(i, :) == 0);
                
                if length(testfi) == 1
                    numlabels(end+1, :) = {label, @(x, h, f, s, ind, i) h(baseindex+testfi, i)};
                    plotlabels(end+1, :)= {label, @(x, h, f, s, ind, i) h(baseindex+testfi, ind)};
                    omap.h{baseindex+testfi} = label;
                    
                else
                    for k = 1:length(testfi)
                        ilabel = sprintf('%s%s', label, num2str(k));
                        numlabels(end+1, :) = {ilabel, @(x, h, f, s, ind, i) h(baseindex+testfi(k), i)};
                        plotlabels(end+1, :)= {ilabel, @(x, h, f, s, ind, i) h(baseindex+testfi(k), ind)};
                        omap.h{baseindex+testfi(k)} = ilabel;
                        
                    end
                    
                end
                
            end
            if ~isempty(numlabels)
                numconfig.declareCategory('tf', 'Test Functions',1e3, false)
                plotsel.declareCategory('tf', 'Test Functions',1e3)
                numconfig.setLabels('tf', numlabels);
                for k = 1:size(plotlabels, 1)
                    plotsel.declareSubCategoryAltName('tf', plotlabels{k, 1}, regexprep(plotlabels{k, 1}, '^test_', '', 'once'), plotlabels{k, 2});
                end
            end
            
        end
        
        
    end
    
end
