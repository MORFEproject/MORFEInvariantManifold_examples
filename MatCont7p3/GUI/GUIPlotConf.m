classdef GUIPlotConf < handle

    
    properties

        configs
        plotselection
        axeshandle = [];
        dimension;
        autofit = false;
        
        windowmanager = [];
        solutionhandler = [];
        pointloader = [];
        lineconfig = [];
        configlauncher = [];
        
        figureposition = [];
    end
    
    properties(Constant)
        AXTITLE = {'Abscissa', 'Ordinate' ,'Applicate'};
    end
    
    events
       selectionChanged 
    end
    
    methods

        function obj = GUIPlotConf(dimension, windowmanager)
            obj.dimension = dimension;
            obj.configs{1} = struct('region', [0 1], 'cat', [], 'subcat', [], 'item', []);
            obj.configs{2} = struct('region', [0 1], 'cat', [], 'subcat', [], 'item', []);
            obj.configs{3} = struct('region', [0 1], 'cat', [], 'subcat', [], 'item', []);
            
            if nargin < 2
                obj.windowmanager = GUIWindowManager([]);
            else
               obj.windowmanager = windowmanager; 
            end
        end
        
        function getConfigLauncher(obj)
           obj.configlauncher; 
        end

        function configurePlotSelection(obj, plotselection)
           obj.plotselection = plotselection; 
           
        end
        
        function setSolutionHandler(obj, sh)
           obj.solutionhandler = sh; 
        end
        function setPointloader(obj, pl)
            obj.pointloader = pl;
        end
        function setLineConfig(obj, lineconfig, configlauncher)
           obj.lineconfig = lineconfig; 
           obj.configlauncher = configlauncher;
        end
        
        function plotSolution(obj)
            if ~isempty(obj.solutionhandler) && ~isempty(obj.solutionhandler.solution)
                oi = obj.solutionhandler.solution.compbranch.getOutputInterpreter();
                %cla(obj.axeshandle); %Clear when (re)drawing solution
                out = obj.getOutputter(oi.getPlotter(), obj.getPointLoaderLabel(obj.solutionhandler.solutionpath));
                out.plotSolution(obj.solutionhandler.solution);
            end
        end
        
        function s = getPointLoaderLabel(obj, path)
            if isempty(obj.pointloader)
                s = '';
            else
               s = obj.pointloader.getLabel(path);
            end
            
            
        end
        
        function plotDiagram(obj)
            if ~isempty(obj.solutionhandler) && ~isempty(obj.solutionhandler.diagram)
                diagram = obj.solutionhandler.diagram;
                solutionnames = diagram.getSolutionNames();
                solutions = diagram.getSolutions();
                curvelabels = diagram.getSolutionLabels(obj.pointloader);
                
                cla(obj.axeshandle);
                for k = 1:length(solutions)
                        oi = solutions{k}.compbranch.getOutputInterpreter();
                        [~, ~, plotsel] = oi.interpret(solutions{k}, solutions{k}.settings, solutions{k}.compbranch);
                        if obj.validatePlotSelection(plotsel)
                            out = obj.getPlotOutput(oi.getPlotter(), plotsel, curvelabels{k});
                            out.plotSolution(solutions{k});
                        else
                           fprintf(2, 'plotDiagram: solution rejected, invalid layout (%s)\n', solutionnames{k}); 
                        end
                end
                
            end
        end
        
        
        function fhandle = generateFigure(obj)
            dimstr = ['Plot' num2str(obj.dimension) 'D'];
            fhandle = figure('NumberTitle' , 'off', 'tag', 'matcont', 'Visible', 'off' , 'Name' , dimstr , 'UserData' , [dimstr ' - '], 'DeleteFcn', @(o,e) onFigureClosing(obj, o), 'KeyPressFcn', @(o,e) SessionOutputKeypress(o, e));            

            
            obj.installToFig(fhandle);
            obj.purge();
            obj.updateAxesLocal();     
            if ~isempty(obj.figureposition)
                set(fhandle, 'Position', obj.figureposition);
                movegui(fhandle, 'onscreen');
            end
            set(fhandle, 'Visible', 'on');
            
        end
        
        function figureClosing(obj, fhandle)
            obj.figureposition = get(fhandle, 'Position');
            
        end
        function installToFig(obj, fig)
            obj.axeshandle = axes('Parent' , fig , 'DeleteFcn' , @(o,e) obj.shutdown() );
            obj.installRegionSync(fig); 
            GUIPlotConf.installMenu(fig, obj)
            
        end
        
        function closeFig(obj)
            if ~isempty(obj.axeshandle) && isvalid(obj.axeshandle)
                delete(get(obj.axeshandle, 'Parent'));
            end
        end

        function d = getDimension(obj)
           d = obj.dimension; 
        end
        
        function [left , right] = getRegion(obj,dim)
            region = obj.configs{dim}.region;
            left = region(1);
            right = region(2);
        end
        
        function s = getData(obj)
            label = '<no label>';
            position = [];
            
            if obj.isValid()
                axisparent = get(obj.axeshandle, 'Parent');
                label = get(axisparent, 'Name');
                position = get(axisparent, 'Position');           
            end
            s = GUIPlotConfLite(obj.getDimension(), obj.configs, label, position);

            % % % % %
        end
        
        function restoreData(obj, data)
            obj.configs = data.configs;
            obj.figureposition = data.figureposition; 
        end
        
        
        function name = getAxLabel(obj, dim)
           if  isempty(obj.configs{dim}.cat)
               name = '<none>';
           else
               name = obj.plotselection.getSelectionLabel(obj.configs{dim}.cat, obj.configs{dim}.subcat, obj.configs{dim}.item);
           end
        end
        
        function updateAxesLocal(obj)
            region = [obj.configs{1}.region  obj.configs{2}.region];
            
            xlbl = obj.getAxLabel(1);

            
            
            ylbl = obj.getAxLabel(2);
            xlabel(obj.axeshandle , xlbl);
            ylabel(obj.axeshandle , ylbl);
            
            labels = [xlbl ',' ylbl];
            
            if (obj.getDimension() > 2)
                zlbl = obj.getAxLabel(3);
                zlabel(obj.axeshandle, zlbl);
                region = [region obj.configs{3}.region];
                labels = [labels ',' zlbl];
            end
            axis(obj.axeshandle , region);
            axisparent = get(obj.axeshandle, 'Parent');
            set(axisparent, 'Name' , [get(axisparent,'UserData') labels]);             
            delete(findall(obj.axeshandle, 'Tag', 'ERR_MSG'))
        end
        
        function setRegion(obj, dim, region)
            obj.configs{dim}.region = region;
            obj.updateAxesLocal();
        end
        
        
         function installRegionSync(obj, figurehandle)
            panh = pan(figurehandle);
            set(panh, 'ActionPostCallback' , @(o,e) obj.syncRegionWithAxes());
            zoomh = zoom(figurehandle);
            set(zoomh, 'ActionPostCallback' , @(o,e) obj.syncRegionWithAxes()); 
         end
         function title = getAxTitle(obj,dim)
             title = obj.AXTITLE{dim};
         end
         
         function syncRegionWithAxes(obj)
             if (~isempty(obj.axeshandle))
                 %fprintf(2, 'SYNC CALL\n')
                 newregion = axis(obj.axeshandle);
                 
                 obj.configs{1}.region = newregion(1:2);
                 obj.configs{2}.region = newregion(3:4);
                 if (obj.getDimension() > 2)
                     obj.configs{3}.region = newregion(5:6);
                 end
                 obj.autofit = false;
                 
                 obj.notify('selectionChanged');
             end
         end
         
         
         function shutdown(obj)
             %fprintf(2, 'on Delete (GUIPlotConf)...\n');
         end
         
         function purge(obj)
             if ~isempty(obj.plotselection)
             for dim = 1:obj.dimension
                 if ~obj.plotselection.isValid(obj.configs{dim}.cat, obj.configs{dim}.subcat, obj.configs{dim}.item)
                     obj.configs{dim}.cat = []; obj.configs{dim}.subcat = []; obj.configs{dim}.item = [];
                 end
             end
             end
         end
         
       
         function setDefaultOnEmpty(obj, session)
             if isempty(obj.plotselection); return; end
             emptySelection = 1;
             for i = 1:obj.dimension
                emptySelection = emptySelection && isempty(obj.configs{i}.cat);
             end
             if emptySelection
             
                % % %
                suggestions = obj.plotselection.getDefaultSelection(session);
                
                for k = 1:min(size(suggestions, 1) , obj.dimension) % FIXME % 3D plot met 2 coordinaten crasht. verkeerde assumpty
                    obj.configs{k}.cat = suggestions{k, 1};
                    obj.configs{k}.subcat = suggestions{k, 2};
                    obj.configs{k}.item = suggestions{k, 3};
                    
                end
             end
         end
         
         function b = validatePlotSelection(obj, plotselection)
            b = 1;
            for dim = 1:obj.dimension
               b = b && plotselection.isValid(obj.configs{dim}.cat, obj.configs{dim}.subcat, obj.configs{dim}.item);
            end
         end
         
         
         function b = hasValidSelection(obj)
             b = 1;
             for k = 1:obj.getDimension()
                 b = b && ~isempty(obj.configs{k}.cat);
             end
         end
         
         function b = isValid(obj)
             obj.purge();
             b = isvalid(obj.axeshandle) && obj.hasValidSelection();
             
             if isvalid(obj.axeshandle); obj.updateAxesLocal(); end
             if ~b
                 obj.displayError('The layout has become invalid');
             end
         end
         
         function b = checkValidSelection(obj)
             b = 1;
             for dim = 1:obj.dimension
                 if ~obj.plotselection.isValid(obj.configs{dim}.cat, obj.configs{dim}.subcat, obj.configs{dim}.item)
                     b = 0;
                 end
             end
             b = b && isvalid(obj.axeshandle) && obj.hasValidSelection();
             delete(findall(obj.axeshandle, 'Tag', 'ERR_MSG'));
             if ~b
                 obj.displayError('The layout has become invalid');
             end             
         end          
         

         
         
         
        function b  = isAlive(obj)
            b = ~isempty(obj.axeshandle) && isvalid(obj.axeshandle);
        end

         
        function displayError(obj, msg)
            if ~isvalid(obj.axeshandle); return; end
            d = axis(obj.axeshandle);
            %cla(obj.axeshandle);
        
            if (obj.getDimension() == 2)
                pos = { (d(1) + d(2))/2 , (d(3) + d(4))/2};
            else
                pos = { (d(1) + d(2))/2 , (d(3) + d(4))/2 , (d(5) + d(6))/2 };
            end
            
            text(pos{:},['\bf ' msg] , 'Parent' , obj.axeshandle  , 'Fontsize' , 16 , 'Color' , 'red' , 'HorizontalAlignment' , 'center' , 'Tag', 'ERR_MSG')
            
        end         
         function [categories, selection] = getCategories(obj, dim)
             categories = obj.plotselection.getCategories();
             selection = obj.configs{dim}.cat;
             
         end
          function [subcategories, selection] = getItems(obj, dim)
             subcategories = obj.plotselection.getItems(obj.configs{dim}.cat, obj.configs{dim}.subcat);
             selection = obj.configs{dim}.item;
          end        
         
          function [list, selection, singleton] = getSubcatItems(obj, dim)
              [list, singleton] = obj.plotselection.getSubcatItems(obj.configs{dim}.cat, obj.configs{dim}.item);
              selection = obj.configs{dim}.subcat;
          end
         
         function setCategory(obj, dim, catname)
             if isempty(catname)
                 obj.configs{dim}.cat = [];
                 obj.configs{dim}.subcat = [];
                 obj.configs{dim}.item = [];
             else
                 
                 subcats = obj.plotselection.getSubCategories(catname);
                 subcat = subcats{1};
                 items = obj.plotselection.getItems(catname, subcat);
                 item = items{1};
                 obj.configs{dim}.cat = catname;
                 obj.configs{dim}.subcat = subcat;
                 obj.configs{dim}.item = item;
                 
             end
             obj.updateAxesLocal();
             obj.notify('selectionChanged');
         end
         function setSelection(obj, dim, subcat, item)
             obj.configs{dim}.subcat = subcat;
             obj.configs{dim}.item = item;
             obj.updateAxesLocal();
             obj.notify('selectionChanged');
         end
         function setAutoFit(obj, bool)
             obj.autofit = bool;
         end
         
         
         
         function o = showLayoutWindow(obj)
             o = [];
             if ~isempty(obj.plotselection)
                obj.purge();
                obj.updateAxesLocal();
                fposition = get(get(obj.axeshandle, 'Parent'), 'Position');
                o  = GUIPlotConfigMainPanel(obj.windowmanager, obj, fposition);
             else
                fprintf(2, 'ERROR: no layout available, please select an initial point type in the menu\n'); 
             end
         end
         
         
         function fitAxis(obj)
             axis(obj.axeshandle , 'tight');
             obj.syncRegionWithAxes();
         end
         
         
         function o = getOutputter(obj, plotter, curvelabel)
            if ~obj.checkValidSelection();
                o = [];
                return;
            end
            o = obj.getPlotOutput(plotter, obj.plotselection, curvelabel);

             
         end
         function o = getPlotOutput(obj, plotter, plotselection, curvelabel)
             
             get = @(dim) plotselection.getSelectFunction(obj.configs{dim}.cat, obj.configs{dim}.subcat, obj.configs{dim}.item);
             
             [fx, pointplotx] = get(1);
             [fy, pointploty] = get(2);
             
             fz = [];
             pointplotz = [];
             if obj.getDimension() > 2
                 [fz, pointplotz] = get(3);
             end
             
             if isempty(obj.lineconfig); obj.lineconfig = GUILineOptions(); end 
             
             
             plotconfig = struct();
             plotconfig(1).curve = obj.lineconfig.getLineSetting(plotselection.getSolutionLabel());
             plotconfig(1).specialpoint = obj.lineconfig.getSingPointSetting();
             plotconfig(1).specialcurve = obj.lineconfig.getSingCurveSetting();
             plotconfig(1).label = obj.lineconfig.getSingLabelSetting();
             plotconfig(1).modification = GUICurveModifications.getMod(plotselection.getSolutionLabel(), obj.lineconfig);
             
             
             o = plotter(obj.axeshandle, fx, fy, fz, [pointplotx pointploty pointplotz], plotconfig, obj.pointloader, curvelabel);
             
         end
         
         
         
    end
        
    methods(Static)
        
        function installMenu(fhandle , plotconf)
            mhandle = uimenu(fhandle, 'Label' , 'MatCont');
            uimenu(mhandle , 'Label' , 'Layout' , 'Callback' , @(o,e)  plotconf.showLayoutWindow());
            
            uimenu(mhandle , 'Label' , 'Redraw Curve' , 'Callback' , @(o,e) plotconf.plotSolution() , 'Separator' , 'on');
            uimenu(mhandle , 'Label' , 'Redraw Diagram' , 'Callback' , @(o,e) plotconf.plotDiagram() );
            uimenu(mhandle , 'Label' , 'Clear' , 'Callback' , @(o,e) cla(plotconf.axeshandle));
            uimenu(mhandle , 'Label' , 'Fit Range' , 'Callback' , @(o,e) plotconf.fitAxis());
            
            
            if ~isempty(plotconf.configlauncher) && ~isempty(plotconf.windowmanager)
                GUIWindowLaunchButton(mhandle , plotconf.windowmanager , 'uimenu', ...
                    plotconf.configlauncher , 'coloropt', 'separator', 'on');
            end
        end

    end
end 

function onFigureClosing(conf, figurehandle)
    if isvalid(conf)
       	conf.figureClosing(figurehandle);
    end
end

