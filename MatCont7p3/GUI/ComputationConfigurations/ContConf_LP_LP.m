classdef ContConf_LP_LP < ContConf_LP
    

    methods
        function obj = ContConf_LP_LP()
            obj = obj@ContConf_LP();
        end
        
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), 'LP'));      
       
        end
        function configureSettings(obj, settings)
            configureSettings@ContConf_LP(obj, settings);
            paramsetting = settings.getSetting('parameters');
            paramsetting.revive(settings, true, true); %active:true, branch:true
        end
        
        
        function p = getPrioritynumber(~)
           p = 1;
       end
    end
    methods(Access=protected)
       function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            paramsetting = settings.getSetting('parameters');
            bp = paramsetting.getBranch();
            [x0,v0] = obj.initFunction(handle,x0,param, activeParams, bp);
        end        
    end
end
