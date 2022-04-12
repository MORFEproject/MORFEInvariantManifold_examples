classdef HTHomCurve < ContCurve
    
    
    
    
    methods
         function obj = HTHomCurve(varargin)
             obj = obj@ContCurve(varargin{:}); 
         end
        
        
        
        function settings = switchAtPoint(obj, index)
            settings = switchAtPoint@ContCurve(obj, index);
            
            s_indices = [obj.s.index];
            ii = find(s_indices == index, 1);
            
            if ~isempty(ii) && strcmpi(strip(obj.s(ii).label), 'hthom')
                
                try
                    settings = HTHom_pointload(obj, settings, index);
                catch ERROR
                    errmsg = ERROR.message;
                    if ~isempty(ERROR.identifier)
                        rethrow(ERROR);
                    end
                    
                    errordlg(errmsg, 'Loading HTHom point');
                    return
                end
                
                    

            end
        end
        
        
    end
    methods(Access=protected)
        function l = getFileLoader(~)
            l = @HTHomCurve.loader;
        end
        
    end
    methods(Static)
        function obj = loader(data)
            obj = HTHomCurve(data.gui.settings, data.gui.compbranch, data.x, data.v, data.s, data.h, data.f, data.globals);
        end
    end    
    
    
end
