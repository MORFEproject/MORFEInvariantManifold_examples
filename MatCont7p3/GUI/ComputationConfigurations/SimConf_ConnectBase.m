classdef SimConf_ConnectBase < CompConf
    
    properties
        simconf;
        type;
        name;
        connectinit;
        pointlabel = 'EP';
    end
    
    
    methods
        
        function b = hasTarget(~) 
            b = 'eps1';
        end        
        
        function s = getLabel(obj)
            s = ['connect_' obj.type '_' obj.simconf.getLabel()];
        end
        function p = getPrioritynumber(obj)
            p = obj.simconf.getPrioritynumber();
        end
        function s = toString(obj)
            s = sprintf('%s',obj.getName());
        end
        function l = getSolutionLabel(obj)
            l = ['Connec' obj.type];
        end
        function s = getName(obj)
            s =  obj.name;
        end
        
        function oi = getOutputInterpreter(~)
            oi = CLSimOutputInterpreter();
        end
        
        
        function obj = SimConf_ConnectBase(name, type, simconf, connectinit)
            obj.simconf = simconf;
            obj.type = type;
            obj.name = name;
            obj.connectinit = connectinit;
        end
        
        function b = isAvailable(obj, settings)
            initialpoint = settings.IP;
            b = ~isempty(settings.system) && strcmp(initialpoint.getILabel(), obj.pointlabel);
            
        end
        
        function configureSettings(obj, settings)
            obj.simconf.configureSettings(settings);
            set = settings.getSetting('eventfunction');
            if ~isempty(set); set.setVisible(0); end %disable option for eventfunction.
            obj.install(settings, 'SParamTestTolerance', settings.getVal('TestTolerance', 1e-05), InputRestrictions.POS, [0, 10000, 1]);
        end
        
        function alist = actions(obj)
            alist(1).label = 'Forward';
            alist(1).function = @(settings, solution) obj.compute(settings, true);
            alist(1).valid = @(~,~) 1;
            alist(2).label = 'Backward';
            alist(2).function = @(settings, solution) obj.compute(settings, false);
            alist(2).valid = @(~,~) 1;
            
            alist(3).label = '-'; alist(3).function = [];
            
            alist(4).label = 'Extend';
            alist(4).function = @(settings, solution) obj.extend(solution, settings);
            alist(4).valid = @(settings, solution) ~isempty(solution);
            
        end
        
        function [solution, errmsg, overwrite] = compute(obj, settings, forward)

            [start_point, data, errmsg] =  obj.connectinit(settings);
            if isempty(start_point)
                solution = []; overwrite = 0;
                return;
            end
            newsettings = settings.copy();
            newsettings.coord.set(start_point);
            
            [solution, errmsg, overwrite] = obj.simconf.compute(newsettings, forward);
            solution.settings = settings;
            solution.connectData = data;
            solution.compbranch = obj;
            
        end
        
        function [solution, errmsg, overwrite] = extend(obj, solution , settings)
            [solution, errmsg, overwrite] = obj.simconf.extend(solution, settings);
            solution.compbranch = obj;
        end
        
        function b = isHidden(obj)
            b = obj.simconf.isHidden();
        end
        function m = getMethodName(obj)
           m = obj.simconf.getMethodName(); 
        end
    end
    
    methods(Access=protected)
        function install(~, settings, name, initval, restrict, catdata)
            settings.addSetting(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingsHelp.getHelp(name) ));
        end
    end
end
