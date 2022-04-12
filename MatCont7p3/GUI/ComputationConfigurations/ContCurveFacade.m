classdef ContCurveFacade 
    
    properties
        x = [];
        v = [];
        globals = [];
    end
    
    
    methods
        function obj = ContCurveFacade(contcurve, index)
            obj.x = FacadeMatrix(contcurve.x, index);
            obj.v = FacadeMatrix(contcurve.v, index);
            obj.globals = contcurve.globals;
        end
        
        function restoreGlobals(obj)
            names = fieldnames(obj.globals);
            for i = 1:length(names)
                eval(sprintf('global %s; %s = obj.globals.%s;', names{i}, names{i}, names{i}));
            end
        end
        
    end       
end
