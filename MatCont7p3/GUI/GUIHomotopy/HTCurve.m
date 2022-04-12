classdef HTCurve < ContCurve
    
    
    properties(Hidden)
        pointloader = [];
        
    end
    
    methods
         function obj = HTCurve(varargin)
             obj = obj@ContCurve(varargin{:}); 
         end
        
        
        
        function settings = switchAtPoint(obj, index)
            settings = switchAtPoint@ContCurve(obj, index);
            
            s_indices = [obj.s.index];
            ii = find(s_indices == index, 1);
            
            if ~isempty(ii)
                
                try
                    settings = obj.pointloader(obj, settings, index);
                catch ERROR
                    errmsg = ERROR.message;
                    if ~isempty(ERROR.identifier)
                        rethrow(ERROR);
                    end
                    
                    errordlg(errmsg, 'Loading Homotopy Point');
                    return
                end
                
                    

            end
        end
        
        
    end
    methods(Access=protected)
        function l = getFileLoader(~)
            assert(false, 'implement this\n');
        end
        
    end
    methods(Static)
        function obj = loader(data)
            assert(false, 'implement this\n');
        end
    end    
    
    
end
