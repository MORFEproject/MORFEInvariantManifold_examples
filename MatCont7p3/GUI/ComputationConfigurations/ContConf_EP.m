classdef ContConf_EP < ContConf
    properties
        initFunction
        nrActive
    end
    
 
    methods
        function obj = ContConf_EP() %curvedefinition, label, nrActive, defaultPointType,initFunction)
            obj.curvedefinition = @equilibrium;
            obj.label = 'EP';
            obj.nrActive = 1;
            obj.defaultPointType = 'EP';
            obj.initFunction = @init_EP_EP;
            obj.testLabels = {'BP', @(s) ContConf.dimensionCheck(s, 1); 'H', @(s) ContConf.dimensionCheck(s, 2); 'LP', @(s) ContConf.dimensionCheck(s, 1)};
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' func2str(obj.initFunction)];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'eds'};
            
        end

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {'EP', 'BP'}));
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            obj.setupInitialPoint(settings);
            obj.setupTestFunctions(settings, obj.testLabels);
            obj.install(settings, 'eigenvalues', true, InputRestrictions.BOOL, [2, 10, 1]);
            
            
        end
        function options = constructContOptions(obj, options, settings)
            options = constructContOptions@ContConf(obj, options, settings);
            options.Eigenvalues = settings.eigenvalues;
            
        end
        
        function p = getPrioritynumber(~)
            p = 1000;
        end
        
        function msg = check(obj, settings)
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            if length(activeParams) ~= obj.nrActive
                msg = sprintf('You have to select exactly %s free parameter', num2str(obj.nrActive));
            else
               msg = ''; 
            end
        end
    end
    
    
    methods(Access=protected)
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            [x0,v0] = obj.initFunction(handle,x0,param, activeParams);
        end
        
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            valid = true;
            msg = '';
            if length(activeParams) ~= obj.nrActive
                valid = false;
                msg = sprintf('You have to select exactly %s free parameter', num2str(obj.nrActive));
            end
        end

        
        function s = postProcess(~, s)  %EP specific
            for i = 1:length(s)
                if strcmp(strip(s(i).label), 'H') && contains(lower(s(i).msg), 'neutral saddle')
                    s(i).internallabel = 'NE';
                end
            end
        end
        

    end
end
