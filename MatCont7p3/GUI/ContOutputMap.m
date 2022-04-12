classdef ContOutputMap
    properties
        x = {}
        %v = {}
        s = {}
        h = []
        f = {}
    end
    
    methods
        function obj = ContOutputMap(xdim, vdim, hdim, fdim)
            obj.x = cell(1, xdim);
            %obj.v = cell(1, vdim);
            obj.h = cell(1, hdim);
            obj.f = cell(1, fdim);
        end
        
        
        function disp(obj)
            fprintf('-OUTPUT-INTERPRET---\n');
            for z = ['x', 'h', 'f']
               for i = 1:length(obj.(z)) 
                    fprintf('%s(%i): %s\n', z, i, obj.(z){i})
                
               end
               fprintf('\n');
            end
        end
    end
    
    
    
    
end