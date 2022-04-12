classdef ContConf_LC_PD < ContConf_LC
    properties
    end
    
 
    methods
        function obj = ContConf_LC_PD() 
            obj = obj@ContConf_LC();
            obj.initFunction = @init_PD_LC2;
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_LC(obj) '_PD'];
        end
        


        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), {'PD'})); 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_LC(obj, settings);
            obj.setupInitialPoint(settings, false);
            obj.install(settings, 'amplitude', 1e-6, InputRestrictions.POS, [2, 5, 1]);
            
            parammodel = settings.getSetting('parameters');
            parammodel.setEditable(false, true); %values: uneditable, active: editable
        end

        
        function p = getPrioritynumber(~)
            p = 1;
        end
    end
    
    
    methods(Access=protected)
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            [x0,v0] = obj.initFunction(handle, IP.source.x, sdata, activeParams, settings.ntst, settings.ncol, settings.amplitude);
        end
        
        
        
    end
    
end
