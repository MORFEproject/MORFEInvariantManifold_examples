classdef ContConf_HSN_LC < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_HSN_LC() 
            obj.curvedefinition = @homoclinicsaddlenode;
            obj.label = 'HSN';
            obj.defaultPointType = 'HSN';
            obj.testLabels = {'NCH',  @(s) ContConf.dimensionCheck(s, 1)};
            obj.initFunction = @init_LC_HSN;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
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

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b = any(strcmp(initialpoint.getILabel(), {'LC',}));
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        
        function configureSettings(obj, settings)
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
            initial = settings.getSetting('IP');
            newT = initial.currentpointdata.data.T / 2;
            
            settings.custom_T.set(newT);
            if settings.custom_eps0 == 0
                settings.custom_eps0.set(0.01);
            end
            if settings.custom_eps1 == 0
                settings.custom_eps1.set(0.01);
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
        
        %function s = postProcess(~, s)  -> dont.
        %    for i = 2:length(s)-1
        %        s(i).internallabel = ['HSN_' strip(s(i).label)];
        %    end
        %end
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
%function [x,v] = init_LC_HSN(odefile, x, s, p, ap, ntst, ncol,extravec,T,eps0,eps1)
            [x0, v0] = obj.initFunction(handle, IP.source.x, sdata, param, activeParams, settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
            
            global homds;
            settings.custom_T.set(homds.T);
            
        end        
    end
    
    
end
