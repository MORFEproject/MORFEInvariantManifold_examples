classdef ContConf_HSN_HTHSN < ContConf_HSN_HSN

 
    methods
        function obj =  ContConf_HSN_HTHSN() 
            obj = obj@ContConf_HSN_HSN();
            obj.initFunction = @init_HTHSN_HSN;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_HSN_HSN(obj) '_HTHSN'];
        end
        

        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            plabel = initialpoint.getILabel();
            b =  any(strcmp(plabel, 'HTHSN'));
        end
        

        function p = getPrioritynumber(~)
            p = 1000000;
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf_HSN_HSN(obj, settings, false);
            
            homoclinicparameters = settings.getSetting('homoclinicparameters');
            htds = settings.htds;
            homoclinicparameters.setValue([htds.T, htds.eps0, htds.eps1]);

            
        end
    end

    methods(Access=protected)     
        function [x0, v0, errmsg] = performInit(obj, handle, point, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            hparam = settings.getSetting('homoclinicparameters');
            extravec = hparam.parametersfree;
            
            x0 = IP.source.x(:, sdata.index);
            v0 = IP.source.v(:, sdata.index);
            x0 = [x0(1:(settings.ntst*settings.ncol+1)*length(point)); point(:)];
            [x0, v0] = obj.initFunction(handle, x0, v0, sdata, param, activeParams, ...
                settings.ntst, settings.ncol, extravec, settings.custom_T, settings.custom_eps0, settings.custom_eps1);
        end
    end
    
    
end
