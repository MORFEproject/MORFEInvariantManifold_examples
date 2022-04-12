classdef ContConf_LP < ContConf_EP
    
    methods
        function obj = ContConf_LP()
            obj = obj@ContConf_EP();
            obj.curvedefinition = @limitpoint;
            obj.label = 'LP';
            obj.nrActive = 2;
            obj.defaultPointType = 'LP';
            obj.initFunction = @init_LP_LP;
            obj.testLabels = {'BT', @(s) ContConf.dimensionCheck(s, 2); 'ZH', @(s) ContConf.dimensionCheck(s, 3); 'CP', @(s) ContConf.dimensionCheck(s, 1)};
        end
        
        
        function b = isAvailable(~, settings)
            %initialpoint = settings.IP;
            %b = any(strcmp(initialpoint.getILabel(), {'BP', 'CP', 'ZH', 'BT'}));      
            b = 0;
        end
        
        function p = getPrioritynumber(~)
           p = 9;
        end
        function list = getGlobalVars(obj)
            list = {'cds', 'lpds'};
            
        end
    end
    
end
