classdef GUIEmptySpace
    
    properties
       width;
       height;
    end
    
    methods
        function obj = GUIEmptySpace(width, height)
           obj.width = width;
           obj.height = height;
        end
        
        function e = Extent(obj)
            e = [0 0 obj.width obj.height];
        end       
        
        function set(~, varargin)
        end
        
    end
    
end