classdef ContConf < CompConf
    %continuation configuration
    properties
        testLabels = {}; %contains testfunctions
        curvedefinition %link to the CL curvedefinitionfile (neimarksacker, etc)
        label %curve label
        defaultPointType %base point type
        contfunction = @cont  %set continuer on default, use this for enabling experimental continuers
        solutionfunction = @ContCurve
    end
    
    methods
        function s = getCurveLabel(obj)
            s = obj.label;
        end
        
        function label = getDefaultPointType(obj)
            label = obj.label;
        end
        
        
        function configureSettings(obj, settings)
            configureSettings@CompConf(obj, settings);

            %Configure settings for continuer:
            obj.install(settings, 'forward', true, InputRestrictions.BOOL, [0, 1, 1]);
            settings.setVisible('forward', false);
            obj.install(settings, 'InitStepsize', DefaultValues.CONTINUERDATA.InitStepSize , InputRestrictions.POS, [1, 1, 1]);
            obj.install(settings, 'MinStepsize', DefaultValues.CONTINUERDATA.MinStepSize , InputRestrictions.POS, [1, 1,2]);
            obj.install(settings, 'MaxStepsize',  DefaultValues.CONTINUERDATA.MaxStepSize, InputRestrictions.POS, [1, 1,3]);
            
            obj.install(settings, 'MaxNewtonIters',  DefaultValues.CONTINUERDATA.MaxNewtonIters, InputRestrictions.INT_g0, [1,2,1]);
            obj.install(settings, 'MaxCorrIters',  DefaultValues.CONTINUERDATA.MaxCorrIters, InputRestrictions.INT_g0, [1,2,2]);
            obj.install(settings, 'MaxTestIters',  DefaultValues.CONTINUERDATA.TestIters, InputRestrictions.INT_g0, [1,2,3]);
            obj.install(settings, 'VarTolerance',  DefaultValues.CONTINUERDATA.VarTolerance, InputRestrictions.POS_0, [1,2,4]);
            obj.install(settings, 'FunTolerance',  DefaultValues.CONTINUERDATA.FunTolerance, InputRestrictions.POS_0, [1,2,5]);
            obj.install(settings, 'TestTolerance',  DefaultValues.CONTINUERDATA.TestTolerance, InputRestrictions.POS_0,  [1,2,6]);
            obj.install(settings, 'Adapt',  DefaultValues.CONTINUERDATA.Adapt, InputRestrictions.INT_ge0, [1,2,7]);
            
            obj.install(settings, 'MaxNumPoints',  DefaultValues.CONTINUERDATA.MaxNumPoints, InputRestrictions.INT_g0, [1,3,1]);
            obj.install(settings, 'CheckClosed',  DefaultValues.CONTINUERDATA.ClosedCurve, InputRestrictions.INT_ge0, [1,3,2]);
            settings.addSetting('option_tsearchorder', CLSetting('TSearchOrder', 1, InputRestrictions.BOOL, 1, 5, 1, CLSettingsHelp.getHelp('option_tsearchorder') ));
            settings.addSetting('option_moorepenrose', CLSetting('Moore-Penrose', 1, InputRestrictions.BOOL, 1, 5, 2, CLSettingsHelp.getHelp('option_moorepenrose') ));
            settings.addSetting('option_increment', CLSetting('Jacobian Increment', 1e-5, InputRestrictions.POS, 1, 5, 3,  CLSettingsHelp.getHelp('option_increment')));
            %Set up user-function settings.
            obj.setupUF(settings);
        end
        
        function options = constructContOptions(obj, options, settings)
            %configure settings to a MATLAB-struct to pass to the continuer
            options.Backward = ~settings.forward;
            
            options.InitStepsize = settings.InitStepsize;
            options.MinStepsize = settings.MinStepsize;
            options.MaxStepsize = settings.MaxStepsize;
            
            options.MaxCorrIters = settings.MaxCorrIters;
            options.MaxNewtonIters = settings.MaxNewtonIters;
            options.MaxTestIters = settings.MaxTestIters;
            
            options.FunTolerance = settings.FunTolerance;
            options.VarTolerance = settings.VarTolerance;
            options.TestTolerance = settings.TestTolerance;
            
            options.Adapt = settings.Adapt;
            options.CheckClosed = settings.CheckClosed;
            options.MaxNumPoints = settings.MaxNumPoints;
            
            options.Increment = settings.option_increment;
            options.MoorePenrose = settings.option_moorepenrose;
            options.TSearchOrder = settings.option_tsearchorder;
            
            %create test-function option
            IgnoreSingularity = [];
            for index = 1:size(obj.testLabels, 1)
                testlabel = obj.testLabels{index, 1};
                internallabel = ['test_', obj.getCurveLabel() , '_', testlabel ];
                if ~settings.getValue(internallabel)
                    IgnoreSingularity(end+1) = index;
                end
            end
            options.IgnoreSingularity = IgnoreSingularity;
            options.Singularities = length(IgnoreSingularity) < length(obj.testLabels); %disable if none selected
            

            %pass along userfunction-settings, if any
            ufsetting = settings.getSetting('userfunctions');
            if ~isempty(ufsetting)
                ufdatamodel = ufsetting.getValue();
                options.UserfunctionsInfo = ufdatamodel.ufdata;
                options.Userfunctions = any([ufdatamodel.ufdata.state]);
            end
            
            
        end

       
        function [labels, internallabels] = getTestLabels(obj)
            %Obtain labels and iternal labels of testfunctions, as placed in settings
            if isempty(obj.testLabels)
                labels = {}; internallabels = {};
               return;
            end
            
            labels = obj.testLabels(:, 1);
             internallabels = cell(1, size(obj.testLabels, 1));
            for index = 1:size(obj.testLabels, 1)
                internallabels{index} = ['test_', obj.getCurveLabel() , '_', labels{index}];
            end
        end
        
        
        function alist = actions(obj)
            %generate the actions for this computaiton

            %each action has a label, a function that tests availability and a function that executes the action
            alist(1).label = 'Forward';
            alist(1).function = @(settings, solution) obj.compute(settings, true);
            alist(1).valid = @(~,~) 1;

            alist(2).label = 'Backward';
            alist(2).function = @(settings, solution) obj.compute(settings, false);
            alist(2).valid = @(~,~) 1;
            
            alist(3).label = '-'; alist(3).function = [];  %add separator
            
            alist(4).label = 'Extend';
            alist(4).function = @(settings, solution) obj.extend(solution);
            alist(4).valid = @(settings, solution) ~isempty(solution);
            
        end        
        function list = getGlobalVars(obj)
            list = {'cds'};
            
        end
        function s = getLabel(~)
            s = 'cont';
        end
        function s = getSolutionLabel(obj)
           s = obj.label; 
        end
        function oi = getOutputInterpreter(~) %Default output interpreter for continuation curves
            oi = CLContOutputInterpreter();
        end
        
        function s = getName(obj)
            if isfield(DefaultValues.CURVENAMES, obj.label)
                s = DefaultValues.CURVENAMES.(obj.label);
            else
                s = obj.label;
            end
        end
        function s = toString(obj)
            s = sprintf('%s (%s)', obj.getName(), obj.getInitStr());
        end
        function s = getInitStr(obj)
            s = '~-~';
        end
        

        %call extend
        function [solution, errormsg, overwrite] = extend(obj, c)
            errormsg = ''; overwrite = 1;  %overwrite previous solution
            c.restoreGlobals();
            [x, v, s, h, f] = obj.contfunction(c.x, c.v, c.s, c.h, c.f, c.globals.cds);
            s = obj.postProcess(s);
            solution = obj.solutionfunction(c.settings, obj, x, v, s, h, f);
        end
        
        %call continuer
        function [solution, errormsg, overwrite] = compute(obj, settings, forward)
            solution = []; overwrite = 0;  %don't overwrite, create new solution
            settings.setValue('forward', forward);
            
            parametersmodel = settings.getSetting('parameters');
            activeParams = parametersmodel.getActive();
            
            system = settings.system;
            x0 = settings.coord;
            x0 = x0(:); %make column vector
            param = settings.parameters; %??
            param = param(:); %make column vector
            
            
            %restore global variables of source curve, if present
            initialPoint = settings.getSetting('IP');
            if ~isempty(initialPoint.source)
               initialPoint.source.restoreGlobals(); 
            end
            
            
            
            %check if settings are valid
            [valid, errormsg] = obj.sanityCheck(system.handle, x0, param, activeParams, settings);
            if ~valid; return; end

            %do init
            [x0, v0, errormsg] = obj.performInit(system.handle, x0, param, activeParams, settings);
            if isempty(x0)
                if isempty(errormsg); errormsg = 'init failed'; end
                return;
            end
            
            %create continuer-struct
            options = contset;
            options.ActiveParams = activeParams;
            options = obj.constructContOptions(options, settings);
            
            %do the actual continuation
            [x, v, s, h, f] = obj.contfunction(obj.curvedefinition, x0, v0, options);
            if isempty(x); errormsg = ''; return; end
            
            %call postProcess (this might add internal labels for use in GUI)
            s = obj.postProcess(s);

            %create solution and return
            solution = obj.solutionfunction(settings, obj, x, v, s, h, f);
            
        end
        
        function bool = hasPeriod(~)
           bool = 0; 
        end
        
    end
    
    methods(Access=protected)
        function install(~, settings, name, initval, restrict, catdata)
            settings.addSetting(name, CLSetting(name, initval, restrict, catdata(1), catdata(2), catdata(3), CLSettingsHelp.getHelp(name) ));
        end
        
        function setupUF(~, settings)
            %place Userfunctions (if any) in settings for selection
            system = settings.getValue('system');
            ufdata = system.getUFData();
            if ~isempty(ufdata)
                ufdatamodel = CLUFdatamodel(ufdata);
                settings.addSetting('userfunctions', ufdatamodel);
                
                for i = 1:length(ufdata)
                    if ufdata(i).valid
                        internallabel = ['uf_' ufdata(i).label];
                        settings.addSetting(internallabel, CLSettingUserFunction(ufdatamodel, i));
                    end
                end
            end
        end
        
        function s = postProcess(~, s)
        end
        
        function setupInitialPoint(~, settings, displaycoords)
            %create 'initial point' in settings.
            if nargin < 3
               displaycoords = true;  %display coordinates of initial point (not relevant for cycles)
            end
            
            system = settings.system;
            %install 'coords:'
            coordinates = CLSettingCoordinates(settings, system.getCoordinates());
            settings.addSetting('coord', coordinates);
            coordinates = settings.getSetting('coord');
            coordinates.setVisible(displaycoords);
            
            parameters = settings.getSetting('parameters');
            if isempty(parameters)
                parameters = CLSettingParameters(settings, system.getParameters(), true); %true: free parameters
                settings.addSetting('parameters', parameters);
            else
                parameters.revive(settings, true); %true: select free parameters (parameters might exist but set on non-free)
            end
            
        end
        
        function setupTestFunctions(obj, settings, labels, defaultvalue)
            if nargin < 4
               defaultvalue = true; 
            end
            
            
            %add testfunctions to settings
            for index = 1:size(labels, 1)
                
                
                testlabel = labels{index, 1};
                validfunction = labels{index, 2}; %check if testfunction can be enabled (dimension too low?)
                internallabel = ['test_', obj.getCurveLabel() , '_', testlabel ];
                
                testlabel_clean = regexprep(testlabel, '^test_', '', 'once');
                if isfield(DefaultValues.TESTFUNCTIONS, testlabel)
                    fullname = DefaultValues.TESTFUNCTIONS.(testlabel);
                else
                    fullname = testlabel_clean;
                end
                
                name = sprintf('%s (%s)', fullname , testlabel_clean);
                setting = CLSetting(name, defaultvalue, InputRestrictions.BOOL, 2, 6, index, CLSettingsHelp.getHelp(internallabel));
                settings.addSetting(internallabel, setting);
                if ~validfunction(settings)
                    setting.setValue(false);
                    setting.setVisible(false);
                end
            end
            
        end
        
        function [x0, v0, errormsg] = performInit(obj, handle, x0, param, activeParams, settings)
            x0= []; v0 = []; errormsg = '';
        end
        function [valid, msg] = sanityCheck(obj, handle, x0, param, activeParams, settings)
            valid = false; msg = 'shutdown';
        end
        
    end
    
    methods(Static)
        function b = dimensionCheck(settings, dim)
            system = settings.system;
            b = system.getDim() >= dim;
        end
    end
end
