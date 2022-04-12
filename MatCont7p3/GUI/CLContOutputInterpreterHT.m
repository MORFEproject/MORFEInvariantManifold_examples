%output interpreter for HTHom
classdef CLContOutputInterpreterHT < CLContOutputInterpreter
    
    properties
       heteroclinic; 
    end
    
    methods
        function obj = CLContOutputInterpreterHT(varargin)
            %set heteroclinic mode if an argument is passed
            obj.heteroclinic = nargin ~= 0;
        end
        
        function [coord, param, period, meta] = getPoint(obj, solution, index)
            meta = struct();
            settings = solution.settings;
            parameters = settings.getSetting('parameters');
            param = parameters.getValue();
            ap = parameters.getActive();
            dim = length(settings.coord);
            [ntst, ncol] = solution.getDiscretization();

            jump = (ncol*ntst + 1)*dim;

            htds = settings.htds;
            if ~isempty(ap)
                coord = solution.x(jump+1:jump+dim,index);
                jump = jump + dim + 1;
                
                if obj.heteroclinic
                    jump = jump + dim; %jump over heteroclinic target point.
                end
                
                for k = 1:length(ap)
                    param(ap(k)) = solution.x(jump, index);
                    jump = jump + 1;
                    
                end
                
                
            else
               coord = htds.x0; 
               jump = jump + 1;
            end
            
            
            jump = jump + length(find(htds.gds.UParamsFree)) + length(find(htds.gds.SParamsFree));
            
            if htds.gds.extravec(1) ~= 0
               period =  solution.x(jump, index);
            else
               period =  htds.T;
            end
            
            
        end
        


        

        function omap = doX(obj, settings, compbranch, coordinates, parameters, ntst, ncol, dim, omap, numconfig, plotsel)
            [index, omap] = doCycle(obj, coordinates, ntst, ncol, dim, omap, numconfig, plotsel);
            
            htds = [];
            if ~isempty(settings.getSetting('htds'))
                htds = settings.htds;
                if ~isfield(htds, 'gds'); return; end            
            else
                return;
            end
            
            %see init_HTxxx_HTxxx, saddle is present if curve not HTHSN and
            %index != 1
            
            hasSaddles = ~contains(compbranch.getCurveLabel(), 'HSN') && compbranch.nextIndex(settings) ~= 1;
            
            
            if ~isempty(parameters.getActive()) || hasSaddles
                numconfig.declareCategory('coordinates_saddle_start', 'Coordinates (saddle)',1, false);
                plotsel.declareCategory('coordinates_saddle_start', 'Coordinates (saddle)', 1);
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
                    numconfig.declareCategory('coordinates_saddle_end', 'Coordinates (target saddle)',1, false);
                    plotsel.declareCategory('coordinates_saddle_end', 'Coordinates (target saddle)', 1);
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
                if ~isempty(parameters.getActive())
                    [index, omap] = obj.doParameters(index, settings, compbranch, parameters, omap, numconfig, plotsel);
                end
            end
            

            curvelabel = compbranch.getCurveLabel();
            if contains(curvelabel, 'Het')
               header = 'Heteroclinic parameters';
            else
               header = 'Homoclinic parameters';
            end
            

            
            UParams = htds.gds.UParams;
            SParams = htds.gds.SParams;
            activeUParams = find(htds.gds.UParamsFree);
            activeSParams = find(htds.gds.SParamsFree);
            if length(activeUParams) + length(activeSParams) > 0  %has active SUParams
                labels = cell(length(activeUParams) + length(activeSParams), 2);
                labelindex = 1;
                
                numconfig.declareCategory('con_params', 'Connection parameters',1e1+4, true);
                plotsel.declareCategory('con_params', 'Connection parameters', 1e1+4);                
                for k = activeUParams
                    paramname = UParams{k, 1};
                    omap.x{index} = paramname;
                    labels(labelindex, :) = {paramname, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('con_params', paramname, @(x, h, f, s, ind, i) x(index, ind));
                    labelindex = labelindex + 1;
                    index = index + 1;
                end
                 for k = activeSParams
                    paramname = SParams{k, 1};
                    omap.x{index} = paramname;
                    labels(labelindex, :) = {paramname, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('con_params', paramname, @(x, h, f, s, ind, i) x(index, ind));
                    labelindex = labelindex + 1;
                    index = index + 1;
                end              
                numconfig.setLabels('con_params', labels);
            end
            
            extravec = htds.gds.extravec;
            activeT = extravec(1); activeEps1 = extravec(3);
            
            if activeT + activeEps1 > 0
                labels = cell(activeT + activeEps1, 2);
                labelindex = 1;    
                
                numconfig.declareCategory('hom_params', header,1e1+8, true);
                plotsel.declareCategory('hom_params', header, 1e1+8);                  
                
                if activeT
                    paramname = 'T';
                    omap.x{index} = paramname;
                    labels(labelindex, :) = {paramname, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('hom_params', paramname, @(x, h, f, s, ind, i) x(index, ind));
                    labelindex = labelindex + 1;
                    index = index + 1;                    
                end
                if activeEps1
                    paramname = 'eps1';
                    omap.x{index} = paramname;
                    labels(labelindex, :) = {paramname, @(x, h, f, s, ind, i) x(index, i)};
                    plotsel.declareSubCategory('hom_params', paramname, @(x, h, f, s, ind, i) x(index, ind));
                    labelindex = labelindex + 1;
                    index = index + 1;                    
                end                
                numconfig.setLabels('hom_params', labels);
            end
            
 
        end
        
        
   
        
        
    end
    
end
