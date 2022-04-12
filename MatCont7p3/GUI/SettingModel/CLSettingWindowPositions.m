classdef CLSettingWindowPositions < CLSettingInterface
   
    properties
        positions;
    end
    
    methods
        
        function obj = CLSettingWindowPositions()
            obj.positions = struct();
            
        end
        
        function setDefault(obj, label, position)
            if ~isfield(obj.positions, label)
               obj.positions.(label) = position; 
            end
            
        end
        function setPosition(obj, label, position)
            obj.positions.(label) = position; 
            
        end
        function pos = getPosition(obj, label)
            pos = obj.positions.(label);
        end
        
        
        function newobj = copy(obj, ~)
            newobj = CLSettingWindowPositions();
            newobj.positions = obj.positions;
        end        
        
    end
end