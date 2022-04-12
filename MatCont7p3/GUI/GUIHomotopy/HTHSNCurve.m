classdef HTHSNCurve < HTCurve
    
    methods
         function obj = HTHSNCurve(varargin)
             obj = obj@HTCurve(varargin{:}); 
             obj.pointloader = @HTHSN_pointload;
         end
        
    end
    methods(Access=protected)
        function l = getFileLoader(~)
            l = @HTHSNCurve.loader;
        end
        
    end
    methods(Static)
        function obj = loader(data)
            obj = HTHSNCurve(data.gui.settings, data.gui.compbranch, data.x, data.v, data.s, data.h, data.f, data.globals);
        end
    end    
    
    
end
