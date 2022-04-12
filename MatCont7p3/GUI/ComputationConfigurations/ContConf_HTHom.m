classdef ContConf_HTHom < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_HTHom() 
            obj.curvedefinition = @homotopysaddle;
            obj.label = 'HTHom';
            obj.defaultPointType = 'HTHom';
            
            obj.testLabels = {'HTHom', @(s) ContConf.dimensionCheck(s, 1)};
                         
            obj.initFunction = @init_HTHom_HTHom;
            
            obj.solutionfunction = @HTHomCurve;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'homds', 'HTHomds'};
            
        end
        function index = nextIndex(~, settings)
           index = HTHom_nextIndex(settings); 
        end
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = any(strcmp(initialpoint.getILabel(), obj.defaultPointType));
        end
        
        function bool = hasPeriod(~)
           bool = 1; 
        end
        
        function configureSettings(obj, settings)
            configureSettings@ContConf(obj, settings);
            obj.setupInitialPoint(settings, false);
            obj.setupTestFunctions(settings, obj.testLabels);  
            s = settings.getSetting('test_HTHom_HTHom'); s.setVisible(false);  %make test invisble.
             
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);
            
            
            settings.addSetting('htds', CLSettingHTDS());  %set htds visible or create if not present.
            
            
            obj.install(settings, 'eps1tol', DefaultValues.STARTERDATA.eps1tol, InputRestrictions.POS, [2, 3, 1000]);
            
            obj.install(settings, 'SParamTestTolerance', settings.getVal('TestTolerance', 1e-05), InputRestrictions.POS, [0, 10000, 1]);
          
            
        end
        

        
        function p = getPrioritynumber(~)
            p = 1;
        end
        
        function msg = check(obj, settings)  %no checks
            msg = '';
        end
        
        
        function oi = getOutputInterpreter(~)
            oi = CLContOutputInterpreterHT();
        end
        
    end

    methods(Access=protected)
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
           msg = obj.check(settings);
           valid = isempty(msg);
        end        
        
        function [x0, v0, errmsg] = performInit(obj, handle, x0, param, activeParams, settings)
            errmsg = '';
            x0 = []; v0 = [];

            try
                [x0, v0] = HTHom_performInit(handle, x0, param, activeParams, settings);
            catch ERROR
                x0 = []; v0 = [];
                errmsg = ERROR.message;
                if ~isempty(ERROR.identifier)
                   rethrow(ERROR); 
                end
                
                
                return
            end
            
            
        end        
    end
    
    
end
