classdef ContConf_Hom_NCH < ContConf_Hom_Hom

 
    methods
        function obj = ContConf_Hom_NCH() 
            obj = obj@ContConf_Hom_Hom();
            obj.initFunction = @init_NCH_Hom;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_Hom_Hom(obj) '_NCH'];
        end
        

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getLabel();
            b =  any(strcmp(plabel, 'NCH'));
        end
        
        function p = getPrioritynumber(~)
            p = 1000000;
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
