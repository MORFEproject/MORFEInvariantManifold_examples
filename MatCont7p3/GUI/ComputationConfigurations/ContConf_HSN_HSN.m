classdef ContConf_HSN_HSN < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_HSN_HSN() 
            obj.curvedefinition = @homoclinicsaddlenode;
            obj.label = 'HSN';
            obj.defaultPointType = 'HSN';
            obj.testLabels = {'NCH',  @(s) ContConf.dimensionCheck(s, 1)};
            obj.initFunction = @init_HSN_HSN;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label '_' obj.defaultPointType];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
         function b = isCycle(obj, settings)
            initial = settings.getSetting('IP');
            b = ~isempty(initial.currentpointdata) && ~isempty(initial.currentpointdata.data) && ~isempty(initial.currentpointdata) && ~isempty(initial.currentpointdata.data.T);
        end       
        
        function list = getGlobalVars(obj)
            list = {'cds', 'homds'};
            
        end

        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {obj.defaultPointType,}));
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        
        function configureSettings(obj, settings , ~)
            configureSettings@ContConf(obj, settings);
            obj.setupInitialPoint(settings, false);
            obj.setupTestFunctions(settings, obj.testLabels);
            
            % % % Homoclinic parameters
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);
            
            
            homoclinicparameters = settings.getSetting('homoclinicparameters');
            if isempty(homoclinicparameters)
                homoclinicparameters = CLSettingCustomParameters(settings, {'T', 'eps0', 'eps1'}, true);
                settings.addSetting('homoclinicparameters', homoclinicparameters);
            else
                homoclinicparameters.revive(settings, true);
            end
            
            obj.install(settings, 'eigenvalues', true, InputRestrictions.BOOL, [2, 10, 1]);
            
            if ~obj.isCycle(settings); return; end
            
            if nargin < 3  %two arguments: load from Hom/Hsn, more: do not load T/eps0/eps1
                initial = settings.getSetting('IP');
                newT = initial.currentpointdata.data.T;
                settings.custom_T.set(newT);
                settings.custom_eps0.set(initial.currentpointdata.data.meta.eps0);
                settings.custom_eps1.set(initial.currentpointdata.data.meta.eps1);
            
            end
                      

        end
        
        
        function options = constructContOptions(obj, options, settings)
            options = constructContOptions@ContConf(obj, options, settings);
            options.Eigenvalues = settings.eigenvalues;
        end
        
        function p = getPrioritynumber(~)
            p = 10000;
        end
        
        function msg = check(obj, settings)
            msg = '';
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            
            if length(activeParams) ~= 2
                msg = 'Two free system parameters are needed';
            end
        end
        
        
        function oi = getOutputInterpreter(~)
            oi = CLContOutputInterpreterConnection();
        end

    end

    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
           msg = obj.check(settings);
           valid = isempty(msg);
        end        
        

        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            %init_HSN_HSN(odefile, x, v, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)`

            [x0, v0] = obj.initFunction(handle, IP.source.x, IP.source.v, sdata, param, activeParams, ...
                settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
            
            global homds;
            settings.custom_T.set(homds.T);
            
        end        
    end
    
    
end
