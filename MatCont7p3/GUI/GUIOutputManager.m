classdef GUIOutputManager < handle
    
    properties
        plotlist;
        session;
        numericpanel;
        
        listener;
        numericmodel;
        plotselection;
        
        
    end
    
    methods
        function obj = GUIOutputManager(session)
            obj.plotlist = {};
            obj.numericpanel = [];
            obj.numericmodel = [];
            obj.plotselection = [];
            
            obj.session = session;
            obj.listener = cell(1, 4);
            obj.listener{1} = session.addlistener('settingChanged', @(o, e) obj.onSettingChanged(session));
            obj.listener{2} = session.addlistener('settingsChanged', @(o, e) obj.onSettingsChanged(session));
            obj.listener{3} = session.addlistener('computationChanged', @(o, e) obj.onSettingChanged(session));
            obj.listener{4} = session.addlistener('diagramChanged', @(o, e) obj.onSettingChanged(session));
            obj.onSettingChanged(session);
        end
        
        
        function destructor(obj)
            delete(obj.listener{:});
        end
        
        function onChangeSystem(obj)
            obj.session.windowmanager.closeLayoutWindows();
            obj.session.windowmanager.closeWindow('monitor');
            obj.session.windowmanager.closeWindow('diagramorganizer');
            obj.session.windowmanager.closeWindow('numeric');
            
            delete(obj.numericpanel);
            for k = 1:length(obj.plotlist)
                obj.plotlist{k}.closeFig();
            end
            
            for k = 1:length(obj.plotlist)
                delete(obj.plotlist{k});
            end
            
            obj.plotlist = {};
            obj.numericpanel = [];
        end
        
        
        function onSettingChanged(obj, session)
            if session.isUnlocked()
                oi = session.computation.getOutputInterpreter();
                if ~isempty(oi)
                    [~, obj.numericmodel, obj.plotselection] = oi.interpret([], session.settings, session.computation);

                    for k = 1:length(obj.plotlist)
                        obj.plotlist{k}.configurePlotSelection(obj.plotselection);
                        obj.plotlist{k}.setSolutionHandler(obj.session.solutionhandle);
                    end

                    obj.session.windowmanager.closeLayoutWindows();
                end
            end
        end
        
        function onSettingsChanged(obj, session)
            obj.numericmodel = obj.session.settings.numericconfig;
            if obj.hasNumeric()
                obj.numericpanel.changeNumericModel(obj.numericmodel);
            end
            obj.onSettingChanged(session);
        end
        
        
        function plotconf = recoverPlot(obj, dimension)
            for k = length(obj.plotlist):-1:1
                % Recover previous plot..
                if obj.plotlist{k}.getDimension() == dimension && ~obj.plotlist{k}.isValid()
                    %fprintf('recover plot from memory\n');
                    plotconf = obj.plotlist{k};
                    return
                    
                end
            end
            previous = obj.session.settings.getSetting('previousplots');
            plotconf = previous.recover(dimension, obj.session.getWindowManager());
            if ~isempty(plotconf)
                %fprintf('recover plot from previous\n');
                return;
            end
            
            
            %fprintf('No recovery: new plot\n');
            %no recover, new plot
            plotconf = GUIPlotConf(dimension, obj.session.getWindowManager());
            
            
            
        end
        
        
        function launch(obj, plotdata)
            plotconf = GUIPlotConf(plotdata.dimension, obj.session.getWindowManager());
            plotconf.restoreData(plotdata);
            obj.activatePlot(plotconf, false);
        end
        
        
        
        function plotconf = createPlot(obj, dimension)
            if isempty(obj.plotselection)
                errordlg('No initial point type selected. Please use "Type>Inital Point" to select a type', 'MatCont')
                
            else
                plotconf = obj.recoverPlot(dimension);
                
                obj.activatePlot(plotconf, true);
                
                if dimension == 3
                   view(30 , 30); %WARNING TODO: experimental new
                end
            end
        end
        
        function activatePlot(obj, plotconf, showLayoutWindow)
            oi = obj.session.computation.getOutputInterpreter();
            if ~isempty(oi)
                plotconf.configurePlotSelection(obj.plotselection);
                plotconf.setSolutionHandler(obj.session.solutionhandle);
                plotconf.setPointloader(obj.session.pointloader);
                plotconf.setLineConfig(obj.session.globalsettings.lineoptions, @(fhandle) GUILineConfigPanel(fhandle , obj.session));
                plotconf.setDefaultOnEmpty(obj.session);
            end
            
            plotconf.generateFigure();
            drawnow;
            if plotconf.dimension == 3
                view(30 , 30); %WARNING TODO: experimental new
            end
            valid = plotconf.isvalid();
            if showLayoutWindow && ~plotconf.hasValidSelection()
                errordlg('Please fill in missing values in the layout window', 'Layout is incomplete');
                plotconf.showLayoutWindow();
            end
            
            obj.plotlist{end+1} = plotconf;
        end
        
        function b = hasNumeric(obj)
            b = ~isempty(obj.numericpanel) && isvalid(obj.numericpanel);
        end
        
        function createNumeric(obj, fighandle)
            
            if obj.hasNumeric()
                delete(obj.numericpanel)
            end
            obj.numericpanel = GUINumericPanel(fighandle, obj.numericmodel);
            
        end
        
        
        function configureSessionOutput(obj)
            global sOutput;
            if obj.hasNumeric()
                outputlist = {obj.numericpanel};
            else
                outputlist = {};
            end
            oi = obj.session.computation.getOutputInterpreter();
            
            [~, ~, obj.plotselection] = oi.interpret([], obj.session.settings, obj.session.computation);
            
            purge = [];
            for k = 1:length(obj.plotlist)
                obj.plotlist{k}.configurePlotSelection(obj.plotselection);
                outputter = obj.plotlist{k}.getOutputter(oi.getPlotter(), obj.session.pointloader.getCurrentLabel());
                if ~isempty(outputter)
                    pp = obj.session.settings.getSetting('previousplots');
                    pp.addPlot(obj.plotlist{k});
                    
                    outputlist{end+1} = outputter;
                else
                    if ~obj.plotlist{k}.isAlive()
                        purge(end+1) = k;
                    end
                end
            end
            obj.plotlist(purge) = []; %TODO FIXME: only purge invisibles?
            
            for item = {'PRCenabled', 'dPRCenabled'}
               setting = obj.session.settings.getSetting(item{1});
               if ~isempty(setting)
                   outputter = setting.getOutputter();
                   if ~isempty(outputter); outputlist{end+1} = outputter; end
               end
            end
            
            
            
            sOutput.outputlist = outputlist;
        end
        
        function list = savePlotConfs(obj)
            list = cell(1, 0);
            for k = 1:length(obj.plotlist)
                dat = obj.plotlist{k}.getData();
                if ~isempty(dat.figureposition)
                    list{end+1} = dat;
                end
            end
            
        end
        
        function loadPlotConfs(obj, list)
            for k = 1:length(list)
                obj.launch(list{k});
            end
        end
        
        
        
    end
    
    
end
