classdef ContConf_EP_point < ContConf_EP
    properties
       pointlabel 
        
    end
    methods
        function obj = ContConf_EP_point(pointlabel, initFunction)
            obj = obj@ContConf_EP();
            obj.initFunction = initFunction;
            obj.pointlabel = pointlabel;
        end
        function s = getLabel(obj)
            s = [getLabel@ContConf_EP(obj) '-' obj.pointlabel];
        end        
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = strcmp(initialpoint.getILabel(), obj.pointlabel);            
        end
        

        function p = getPrioritynumber(~)
           p = 100;
       end
    end
    
end
