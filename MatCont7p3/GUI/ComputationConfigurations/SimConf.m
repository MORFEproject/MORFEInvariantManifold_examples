classdef SimConf < CompConf
    %continuation configuration
    
    properties
        methodname
        defaultrefine
        priority
        allowEvents
    end
    
    methods
        function s = getLabel(obj)
            s = ['sim_' obj.methodname];
        end
        function p = getPrioritynumber(obj)
            p = obj.priority;
        end
        function s = toString(obj)
            s = sprintf('%s', obj.getName());
        end
        function s = toStringCL(obj) 
            s = sprintf('%s (%s)',obj.getName(), obj.methodname);
        end
        function l = getSolutionLabel(~)
            l = 'O';
        end
        function s = getName(~)
            s =  'Orbit';
        end
        
        function b = hasTarget(~) 
            b = [];
        end
        function oi = getOutputInterpreter(~)
            oi = CLSimOutputInterpreter();
        end
        function n = getMethodName(obj)
            n = obj.methodname;
        end
        
        function obj = SimConf(methodname, defaultrefine, allowEvents, priority)
            obj.methodname = methodname;
            obj.priority = priority;
            if nargin < 2
                obj.defaultrefine = 1;
            else
                obj.defaultrefine = defaultrefine;
            end
            
            obj.allowEvents = allowEvents;
            
        end
        function b = isAvailable(~, settings)
            initialpoint = settings.IP;
            b = ~isempty(settings.system) && strcmp(initialpoint.getILabel(), 'P');
            
        end
        
        function b = isHidden(obj)
           b = 1; 
        end
        
        function configureSettings(obj, settings)
            configureSettings@CompConf(obj, settings);
            
            obj.install(settings, 'forward', true, InputRestrictions.BOOL, [0, 1, 1]);
            settings.setVisible('forward', false);
            %odemethods = {'ode45', 'ode23', 'ode113', 'ode15s', 'ode23s', 'ode23t', 'ode23tb', 'ode78', 'ode87'};
            %obj.install(settings, 'TMP_Method_TMP', obj.methodname, InputRestrictionCategory(odemethods), [1, 1, 1]);
            
            system = settings.system;
            
            settings.addSetting('time', CLSetting(system.getTimeName(), 0, InputRestrictions.NUM, 2, 1, 0, CLSettingsHelp.getHelp('time') ));
            coordinates = CLSettingCoordinates(settings, system.getCoordinates());
            settings.addSetting('coord', coordinates);
            
            parameters = settings.getSetting('parameters');
            if isempty(parameters)
                parameters = CLSettingParameters(settings, system.getParameters());
                settings.addSetting('parameters', parameters);
            else
                parameters.revive(settings, false); %false: no free parameters 
            end
            
            obj.install(settings, 'Interval', 1, InputRestrictions.NUM, [3, 1, 1]);
            
            settings.addSetting('InitStepSize_sim', CLSettingBlank('InitStepSize', [], InputRestrictions.POS, 3, 1, 3, CLSettingsHelp.getHelp('InitStepSize_sim'), '<automatic>'));
            settings.addSetting('MaxStepSize_sim', CLSettingBlank('MaxStepSize', [], InputRestrictions.POS, 3, 1, 4, CLSettingsHelp.getHelp('MaxStepSize_sim'), '<automatic>'));
            
            
            obj.install(settings, 'RelTolerance', 1e-3, InputRestrictions.POS, [3, 1, 5]);
            obj.install(settings, 'AbsTolerance', 1e-6, InputRestrictions.POS, [3, 1, 6]);  %can also be vector (component wise) FIX?
            obj.install(settings, 'Refine', obj.defaultrefine, InputRestrictions.INT_g0, [3, 1, 7]);  %remark: ODE45 is '4' FIX?
            obj.install(settings, 'Normcontrol', false, InputRestrictions.BOOL, [3, 1, 8]);
            
            settings.addSetting('eventfunction' , CLSettingEventFunction('EventFunction'));
            
            ef = settings.getSetting('eventfunction');
            ef.setVisible(obj.allowEvents);
            
            refinesetting = settings.getSetting('Refine');
            %if ~refinesetting.isAdjusted()
            refinesetting.forceValue(obj.defaultrefine);
            %end
            
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
            solution = []; errmsg = []; overwrite = 0;
            
            settings.setValue('forward', forward);
            
            [tspan, optionsODE] = obj.buildOptions(settings);
            system = settings.system;
            handles = system.handle();
            x0 = settings.coord;
            x0 = x0(:); %make column vector
            param = num2cell(settings.parameters);
            method = str2func(obj.methodname);
            
  
            if isempty( optionsODE.Events )
                [t, y] = method(handles{2}, tspan, x0, optionsODE, param{:});
                tE = []; yE = []; iE = [];
            else
                global Matcont_call
                Matcont_call = @(tspan, x0, opt) method(handles{2}, tspan, x0, opt, param{:});  
                [t, y, tE, yE, iE] = method(handles{2}, tspan, x0, optionsODE, param{:});
            end
            solution = SimCompSolution(settings, obj, t, y, tE, yE, iE,  method, tspan, x0', optionsODE, cell2mat(param));
            
        end
        
        function [solution, errmsg, overwrite] = extend(obj, solution , settings)
            errmsg = []; overwrite = 1;
        %    interval = settings.Interval;
            if settings.forward
             interval = settings.Interval;
            else
            interval = -settings.Interval;
            end
            tspan = [solution.tspan(2), solution.tspan(2) + interval]; %interval is negative when 'backwards' was set.
            x0 = solution.y(end,:)';
            param = num2cell(solution.param);
            method = solution.method;
            optionsODE = solution.options;
            system = solution.settings.system;
            handles = system.handle();
            
            if isempty( optionsODE.Events )
                [t, y] = method(handles{2}, tspan, x0, optionsODE, param{:});
                tE = []; yE = []; iE = [];
            else
                global Matcont_call
                Matcont_call = @(tspan, x0, opt) method(handles{2}, tspan, x0, opt, param{:});                
                [t, y, tE, yE, iE] = method(handles{2}, tspan, x0, optionsODE, param{:});
            end
            
            t = [solution.t; t];
            y = [solution.y; y];
            tE = [solution.tE; tE];
            yE = [solution.yE; yE];
            iE = [solution.iE; iE];
            solution = SimCompSolution(solution.settings, obj, t, y, tE, yE, iE, method, [solution.tspan(1), tspan(2)], x0', optionsODE, cell2mat(param));
        end
        
    end
    
    methods(Access=protected)
        function [tspan, options] = buildOptions(~, settings)
            system = settings.system;
            handles = system.handle();
            
            options = odeset(); %empty
            
            
            options = odeset(options, 'MaxStep', settings.MaxStepSize_sim);   %value is [] when not set in GUI
            options = odeset(options, 'InitialStep', settings.InitStepSize_sim);   %value is [] when not set in GUI
            options = odeset(options, 'RelTol', settings.RelTolerance);
            options = odeset(options, 'AbsTol', settings.AbsTolerance);
            options = odeset(options, 'Refine', settings.Refine);
            options = odeset(options, 'NormControl', CLbool2text(settings.Normcontrol));
            
            options = odeset(options, 'Jacobian', handles{3});
            options = odeset(options, 'Hessians', handles{5});
            
            ef = settings.getSetting('eventfunction');
            if ~isempty(ef) && ef.isVisible() && ~isempty(ef.getValue())
                eventfunction = ef.getValue();
                options = odeset(options, 'Events', eventfunction);
            end
            
            t = settings.time;
            interval = settings.Interval;
            tspan = [t, t+interval];
            if ~settings.forward
                interval = abs(diff(tspan));
                tspan = [tspan(1), tspan(1) - interval];
            end
            
            global sOutput
            if ~isempty(sOutput) && sOutput.enabled
                options = odeset(options, 'OutputFcn', @GUIODEOutput);
                if ~isempty(ef)
                    global Matcont_options
                    Matcont_options = options;
                end
            end
 
        end
        
        
        function install(~, settings, name, initval, restrict, catdata)
            settings.addSetting(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingsHelp.getHelp(name) ));
        end
    end
    
    
end
