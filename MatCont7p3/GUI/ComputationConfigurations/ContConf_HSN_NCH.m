classdef ContConf_HSN_NCH < ContConf_HSN_HSN
    
    
    methods
        function obj = ContConf_HSN_NCH()
           obj =  obj@ContConf_HSN_HSN();
           obj.defaultPointType = 'NCH';
           obj.initFunction = @init_NCH_HSN;
        end
        
    end
    
    
    
end