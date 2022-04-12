classdef ContConf_LP_point < ContConf_LP
    
    properties
       pointlabel 
    end
    methods
        function obj = ContConf_LP_point(pointlabel, initFunction)
            obj = obj@ContConf_LP();
            obj.pointlabel = pointlabel;
            obj.initFunction = initFunction;
        end
        
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel));      
       
        end
        
        function p = getPrioritynumber(~)
           p = 5;
       end
    end
    
end
