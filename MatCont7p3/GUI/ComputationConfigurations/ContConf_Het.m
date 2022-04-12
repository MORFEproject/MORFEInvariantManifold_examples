classdef ContConf_Het < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_Het(pointtype, initfunction) 
            obj.curvedefinition = @heteroclinic;
            obj.label = 'Het';
            obj.defaultPointType = pointtype;
            obj.testLabels = {};         
            obj.initFunction = initfunction;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label '_' obj.defaultPointType];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'hetds'};
            
        end

        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  strcmp(plabel, obj.defaultPointType);
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        

        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            obj.setupInitialPoint(settings, false);
            obj.setupTestFunctions(settings, obj.testLabels, false);
            
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
            
          
            
            if strcmp(obj.defaultPointType, 'HTHet')
                homoclinicparameters = settings.getSetting('homoclinicparameters');
                htds = settings.htds;
                homoclinicparameters.setValue([htds.T, htds.eps0, htds.eps1]);
            elseif strcmp(obj.defaultPointType, 'Het')
                
                initial = settings.getSetting('IP');
                newT = initial.currentpointdata.data.T;
                settings.custom_T.set(newT);
                settings.custom_eps0.set(initial.currentpointdata.data.meta.eps0);
                settings.custom_eps1.set(initial.currentpointdata.data.meta.eps1);
                settings.coord_target.set(initial.currentpointdata.data.meta.coord_target);
            end
        end
        
        
        function options = constructContOptions(obj, options, settings)
            options = constructContOptions@ContConf(obj, options, settings);
            options.Eigenvalues = settings.eigenvalues;
        end
        
        function p = getPrioritynumber(~)
            p = 9999;
        end
        
        function msg = check(obj, settings)
            msg = '';
        end
        
        
        function oi = getOutputInterpreter(~)
            oi = CLContOutputInterpreterConnection('heteroclinic');
        end
        
    end

    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
           msg = obj.check(settings);
           valid = isempty(msg);
        end        
        

        
        function [x0, v0, errmsg] = performInit(obj, handle, ~, param, activeParams, settings)
            
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            x0 = settings.coord;
            x1 = settings.coord_target;
            x0 = [IP.source.x(1:(settings.ntst*settings.ncol+1)*length(x0));x0(:);x1(:)]; 
            v0 = IP.source.v(:,sdata.index);
            [x0, v0] = obj.initFunction(handle, x0, v0, sdata, param, activeParams, settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
            
            global hetds;
            settings.custom_T.set(hetds.T);
            
        end        
    end
    
    
end
