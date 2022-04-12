classdef ContConf_LC_BPC < ContConf_LC
    properties
    end
    
 
    methods
        function obj = ContConf_LC_BPC() 
            obj = obj@ContConf_LC();
            obj.initFunction = @init_BPC_LC;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_LC(obj) '_BPC'];
        end
        


        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {'BPC'})); 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LC(obj, settings);
            

            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(false, false); %values: uneditable, active: uneditable
  
            obj.install(settings, 'amplitude', 1e-6, InputRestrictions.POS, [2, 5, 1]);
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
           
            [x0,v0] = obj.initFunction(handle, IP.source.x, IP.source.v, sdata, settings.ntst, settings.ncol, settings.amplitude);
            
        end
        
        
        
    end
    
end
