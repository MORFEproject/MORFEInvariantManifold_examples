classdef InitialType
%Initial Point/Orbit    
    properties
        publiclabel;
        internallabel;
        name;
        
    end
    methods
        function obj = InitialType(name, publiclabel, internallabel)
            obj.publiclabel = publiclabel;
            obj.name = name;
            if nargin < 3
                obj.internallabel = publiclabel;
            else
                obj.internallabel = internallabel;
            end
           
            
        end
        function label = getLabel(obj)
           label = obj.publiclabel;
        end
        function label = getILabel(obj)
           label = obj.internallabel; 
        end
        function name = getName(obj)
           name = obj.name; 
        end
    end
    
    
    
end
