classdef ContConf_LC_LC < ContConf_LC
    properties
    end
    
 
    methods
        function obj = ContConf_LC_LC() 
            obj = obj@ContConf_LC();
            obj.initFunction = @init_LC_LC;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_LC(obj) '_LC'];
        end
        


        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {'LC', 'BPC', 'NC', 'LPC', 'NS'})); 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LC(obj, settings);
            obj.setupInitialPoint(settings, false);
        end

        
        function p = getPrioritynumber(~)
            p = 10;
        end
    end
    
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
        
            

            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            
            [x0,v0] = obj.initFunction(handle,IP.source.x, IP.source.v, sdata, param, activeParams, settings.ntst, settings.ncol);
            
        end
        
        
        
    end
    
end
