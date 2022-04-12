classdef ContConf_Hom_LC < ContConf_Hom

 
    methods
        function obj = ContConf_Hom_LC() 
            obj = obj@ContConf_Hom();
            obj.initFunction = @init_LC_Hom;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_Hom(obj) '_LC'];
        end
        

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  any(strcmp(plabel, 'LC'));
        end
        
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_Hom(obj, settings);
            
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
        
        
        
        function p = getPrioritynumber(~)
            p = 1000000;
        end
        
        function b = isCycle(obj, settings)
            initial = settings.getSetting('IP');
            b = ~isempty(initial.currentpointdata) && ~isempty(initial.currentpointdata.data) && ~isempty(initial.currentpointdata) && ~isempty(initial.currentpointdata.data.T);
        end
        
        function msg = check(obj, settings)
            if ~obj.isCycle(settings)
               msg = 'Initial Point is not a cycle';
               return;
            end
            
            msg = check@ContConf_Hom(obj, settings);
        end
        
        
    end

    methods(Access=protected)     
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            [x0, v0] = obj.initFunction(handle, IP.source.x, sdata, param, activeParams, settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
        end
    end
    
    
end
