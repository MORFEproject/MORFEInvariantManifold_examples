classdef ContConf_EP_BP < ContConf_EP
    
    methods
        function obj = ContConf_EP_BP()
            obj = obj@ContConf_EP();
            obj.initFunction = @init_BP_EP;
        end
        
        
        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = strcmp(initialpoint.getILabel(), 'BP');            
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf_EP(obj, settings);
            obj.install(settings, 'amplitude', 1e-6, InputRestrictions.POS, [2, 5, 1]);
        end
        function p = getPrioritynumber(~)
           p = 1;
       end
    end
    
    methods(Access=protected)
       function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            IP = settings.getSetting('IP');
            sdata = IP.currentpointdata;
            if isempty(sdata)
               sdata.data = []; %makes init use EP_EP when pointdata is empty. 
            end
            errmsg = '';
            [x0,v0] = obj.initFunction(handle,x0,param, sdata, settings.amplitude);
        end        
    end
end
