classdef ContConf_LC_H < ContConf_LC
    properties
    end
    
 
    methods
        function obj = ContConf_LC_H() 
            obj = obj@ContConf_LC();
            obj.initFunction = @init_H_LC;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_LC(obj) '_H'];
        end
        


        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {'H', 'GH', 'ZH', 'BT', 'HH'})); 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LC(obj, settings);
            obj.setupInitialPoint(settings, true);
            coord = settings.getSetting('coord');
            coord.setVisible(1);
            obj.install(settings, 'amplitude', 1e-6, InputRestrictions.POS, [2, 5, 1]);
        end

        
        function p = getPrioritynumber(~)
            p = 1;
        end
    end
    
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            [x0,v0] = obj.initFunction(handle,x0,param,activeParams, settings.amplitude, settings.ntst, settings.ncol);
            if isempty(x0)
               errmsg = 'Initial point is a Neutral-Saddle'; 
            end
        end
        
        

        
    end
    
end
