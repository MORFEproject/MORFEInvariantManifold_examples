classdef CLCompSwitches < handle
    properties
       switches = {};
    end
    
    
    methods
        function addSwitch(obj, label, description, switcher, data)
            obj.switches{end+1} = struct('label', label, 'description', description, 'switcher', switcher, 'data', data);
        end
        
        function l = length(obj)
            l = length(obj.switches);
        end
        
        function disp(obj)
           for i = 1:length(obj.switches)
                fprintf('%s\n', obj.toString(i));
           end
        end
        
        function s = toString(obj, i)
            s = sprintf('%s: %s (%s)', pad(obj.switches{i}.label,4), pad(obj.switches{i}.description, 40), obj.switches{i}.data.toString(obj.switches{i}.data)); 
        end
        
        function s = getName(obj, i)
            s = sprintf('%s %s', pad(obj.switches{i}.label, 5), obj.switches{i}.description); 
        end
        
        function settings = activateSwitch(obj, index)
           settings = obj.switches{index}.switcher(); 
        end
        
        function data = getData(obj, index)
           data = obj.switches{index}.data ;
            
        end
        
    end
    
    
end