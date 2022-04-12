classdef ContConf_Hom_HTHom < ContConf_Hom

 
    methods
        function obj = ContConf_Hom_HTHom() 
            obj = obj@ContConf_Hom();
            obj.initFunction = @init_HTHom_Hom;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_Hom(obj) '_HTHom'];
        end
        

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  any(strcmp(plabel, 'HTHom'));
        end
        

        function p = getPrioritynumber(~)
            p = 1000000;
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf_Hom(obj, settings);
            
            homoclinicparameters = settings.getSetting('homoclinicparameters');
            htds = settings.htds;
            homoclinicparameters.setValue([htds.T, htds.eps0, htds.eps1]);

            
        end
    end

    methods(Access=protected)     
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            x0 = IP.source.x(:, sdata.index);
            v0 = IP.source.v(:, sdata.index);
            
            [x0, v0] = obj.initFunction(handle, x0, v0, sdata, param, activeParams, settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
        end
    end
    
    
end
