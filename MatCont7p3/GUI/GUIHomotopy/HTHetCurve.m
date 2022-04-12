classdef HTHetCurve < HTCurve
    
    methods
         function obj = HTHetCurve(varargin)
             obj = obj@HTCurve(varargin{:}); 
             obj.pointloader = @HTHet_pointload;
         end
        
    end
    methods(Access=protected)
        function l = getFileLoader(~)
            l = @HTHetCurve.loader;
        end
        
    end
    methods(Static)
        function obj = loader(data)
            obj = HTHetCurve(data.gui.settings, data.gui.compbranch, data.x, data.v, data.s, data.h, data.f, data.globals);
        end
    end    
    
    
end
