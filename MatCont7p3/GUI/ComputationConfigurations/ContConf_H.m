classdef ContConf_H < ContConf_EP
    
    properties
       pointlabel 
    end
    methods
        function obj = ContConf_H(pointlabel, initFunction)
            obj = obj@ContConf_EP();
            obj.curvedefinition = @hopf;
            obj.label = 'H';
            obj.nrActive = 2;
            obj.defaultPointType = 'H';
            
            obj.testLabels = {'BT', @(s) ContConf.dimensionCheck(s, 2); 'ZH', @(s) ContConf.dimensionCheck(s, 3); 'HH', @(s) ContConf.dimensionCheck(s, 4); 'GH', @(s) ContConf.dimensionCheck(s, 2)};
            
            if nargin == 0
                obj.initFunction = @init_H_H;
                obj.pointlabel = 'H';
            else
                obj.initFunction = initFunction;
                obj.pointlabel = pointlabel;     
            end
   
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf_EP(obj) '-' obj.pointlabel];
        end           
        
        function b = isAvailable(obj, settings)
            %system = settings.system; %TODO FIXME
            initialpoint = settings.IP;
            %b = system.getDim() >= 2 && any(strcmp(initialpoint.getILabel(), obj.pointlabel));      
            b = any(strcmp(initialpoint.getILabel(), obj.pointlabel)); 
        end
        
        function p = getPrioritynumber(~)
           p = 6;
        end
       
        function list = getGlobalVars(obj)
            list = {'cds', 'hds'};
            
        end
    end
    
end
