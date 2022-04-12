classdef ContCurve < CompSolution
    
    properties
        x = [];
        v = [];
        s = [];
        h = [];
        f = [];
    end
    
    properties(Hidden)
        globals = struct();
        outputinterpret = [];
        initUsed = '';
    end
    
    methods
        function obj = ContCurve(settings, compbranch, x, v, s, h, f, globals)
            obj = obj@CompSolution(settings, compbranch);
            
            obj.x = x;
            obj.v = v;
            obj.s = s;
            obj.h = h;
            obj.f = f;
            
            obj.initUsed = compbranch.getLabel();
            
            if nargin < 8
                for globalvar = compbranch.getGlobalVars()
                    eval(sprintf('global %s; obj.globals.%s = %s;', globalvar{1}, globalvar{1}, globalvar{1}));
                    eval(sprintf('global %s; %s = []; clear %s', globalvar{1}, globalvar{1}, globalvar{1}));  %TMP, remove global vars as a test
                end
            else
                obj.globals = globals;
            end
            obj.outputinterpret = compbranch.getOutputInterpreter();
            
        end
        
        function restoreGlobals(obj)
            names = fieldnames(obj.globals);
            for i = 1:length(names)
                eval(sprintf('global %s; %s = obj.globals.%s;', names{i}, names{i}, names{i}));
            end
        end
        
        function switches = listSwitches(obj, verbose)
            if nargin < 2
                verbose = false;
            end
            if verbose
               switches = obj.listAllSwitches(); 
            else
                switches = obj.listSpecialSwitches();
            end
        end
        
        function switches = listAllSwitches(obj)    
            switches = CLCompSwitches();
            
            for index = 1:size(obj.x,2)
                [coord, param, period] = obj.outputinterpret.getPoint(obj, index);
                data = struct('npoint', index, 'coord', coord, 'param', param, 'period', period, 'detect', [], 'toString', @(s) [mat2str(s.coord) ' | ' mat2str(s.param)]);
                switches.addSwitch('..', '', @() obj.switchAtPoint(index), data);
            end
        end
        
        function switches = listSpecialSwitches(obj)
            switches = CLCompSwitches();
            for k = 1:length(obj.s)
                [coord, param, period] = obj.outputinterpret.getPoint(obj, obj.s(k).index);
                data = struct('npoint', obj.s(k).index, 'coord', coord, 'param', param, 'period', period, 'detect', obj.s(k), 'toString', @(s) [mat2str(s.coord) ' | ' mat2str(s.param)]);
                switches.addSwitch(strip(obj.s(k).label), obj.s(k).msg, @() obj.switchAtPoint(obj.s(k).index), data);
            end
        end
        
        
        
        function [a, b, c] = interpret(obj)
            [a, b, c] = obj.outputinterpret.interpret(obj, obj.settings, obj.compbranch);
        end
        function m = getStructInterpreter(obj)
            [m, ~, ~] = obj.outputinterpret.interpret(obj, obj.settings, obj.compbranch);
        end
        
        function [b, msg] = showTable(obj, windowname, session)
            if nargin < 3
                session = [];
            end
            b = 1; msg = '';
            if nargin < 2
                windowname = '';
            end
            GUIContCurveTable.previewPanel(obj, windowname, session)
            
            
        end
        
        
        function [ntst, ncol] = getDiscretization(obj)
            ntst = 0; ncol = 0;
            
            names = fieldnames(obj.globals);
            
            for i = 1:length(names)
                if isfield(obj.globals.(names{i}), 'ntst')
                    ntst = obj.globals.(names{i}).ntst;
                    ncol = obj.globals.(names{i}).ncol;
                    return;
                end
            end
        end
        
        function settings = switchAtPoint(obj, index)
            
            s_indices = [obj.s.index];
            
            ii = find(s_indices == index, 1);
            
            settings = obj.settings.copy();
            [coord, param, ~, meta] = obj.outputinterpret.getPoint(obj, index);
            settings.setValue('coord', coord);
            settings.setValue('parameters', param);
            
            default_pt = obj.compbranch.getDefaultPointType();
            
            if ~isempty(ii)
                obj.s(ii).data.meta = meta;
                settings.setValue('IP', obj.s(ii));
                odesys = obj.settings.system;
                if any(strcmp(odesys.getUFLabels(), strip(obj.s(ii).label)))
                    settings.getSetting('IP').assignInternalLabel(default_pt);
                elseif any(strcmp({'99', '00'}, strip(obj.s(ii).label)))
                    settings.getSetting('IP').assignInternalLabel(default_pt);
                    settings.getSetting('IP').assignLabel(default_pt);
                end                
            else
                settings.IP.set(default_pt);
            end
            
            
            settings.getSetting('IP').setSource(ContCurveFacade(obj, index));
            
        end
        
        function save(obj, filename)
            loaderfunction = obj.getFileLoader();
            data = struct('x', obj.x, 'v', obj.v, 's', obj.s, 'h', obj.h, 'f', obj.f, 'globals', obj.globals, 'initUsed', obj.initUsed, 'gui', struct('settings', obj.settings, 'compbranch', obj.compbranch, 'loader', @(dat) loaderfunction(dat)));
            save(filename, '-struct', 'data')
        end
        
        function pr = getPreviewRenderer(obj)
            pr = GUIPreviewRendererContCurve();
        end
        
        

    end
    
    methods(Access=protected)
        function l = getFileLoader(~)
            l = @ContCurve.loader;
        end
        
    end
    methods(Static)
        function obj = loader(data)
            obj = ContCurve(data.gui.settings, data.gui.compbranch, data.x, data.v, data.s, data.h, data.f, data.globals);
        end
    end
    
end
