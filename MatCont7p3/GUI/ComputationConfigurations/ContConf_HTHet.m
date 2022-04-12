classdef ContConf_HTHet < ContConf
    properties
        initFunction;
    end
    
 
    methods
        function obj = ContConf_HTHet() 
            obj.curvedefinition = @homotopyHet;
            obj.label = 'HTHet';
            obj.defaultPointType = 'HTHet';
            obj.testLabels = {'HTHet', @(s) ContConf.dimensionCheck(s, 1)};
            obj.initFunction = @init_HTHet_HTHet;
            
            obj.solutionfunction = @HTHetCurve;
            
        end
        
        function s = getLabel(obj)
            s = [getLabel@ContConf(obj) '_' obj.label];
        end
        
        function s = getInitStr(obj)
            s = func2str(obj.initFunction);
        end
        
        
        function list = getGlobalVars(obj)
            list = {'cds', 'hetds', 'HTHetds'};
            
        end
        function index = nextIndex(~, settings)
           index = HTHet_nextIndex(settings); 
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
            s = settings.getSetting('test_HTHet_HTHet'); s.setVisible(false);  %make test invisble.
             
            obj.install(settings, 'ntst', DefaultValues.STARTERDATA.ntst, InputRestrictions.INT_g0, [2, 4, 1]);
            obj.install(settings, 'ncol', DefaultValues.STARTERDATA.ncol, InputRestrictions.INT_g0, [2, 4, 2]);
            
            settings.addSetting('htds', CLSettingHTDS());  %set htds visible or create if not present.
            obj.install(settings, 'eps1tol', DefaultValues.STARTERDATA.eps1tol, InputRestrictions.POS, [2, 3, 1000]);
            
            obj.install(settings, 'SParamTestTolerance', settings.getVal('TestTolerance', 1e-05), InputRestrictions.POS, [0, 10000, 1]);

            
            
        end

%       function options = constructContOptions(obj, options, settings)
%           options = constructContOptions@ContConf(obj, options, settings);
%           options.Eigenvalues = 1;
%       end
        
        function p = getPrioritynumber(~)
            p = 1;
        end
        
        function msg = check(obj, settings)  %no checks
            msg = '';
        end
        
        
        function oi = getOutputInterpreter(~)
            oi = CLContOutputInterpreterHT('heteroclinic');
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
                [x0, v0] = HTHet_performInit(handle, x0, param, activeParams, settings);
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
