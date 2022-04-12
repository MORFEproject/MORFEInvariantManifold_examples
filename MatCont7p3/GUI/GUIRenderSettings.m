classdef GUIRenderSettings < handle
    properties
        options =  {'Units', 'Pixels', 'Visible', 'on'};
        suggestions = struct('labelsize', 150 );
        
        figurehandle;
        panel;
        session = [];
        
        listener
        idfilter
        
        namemap;
        
        settinglistener = [];
        locklistener;
        enabled;
        windowtitle;
        
        slider;
    end
    
    
    
    methods
        function obj = GUIRenderSettings(fhandle, settings, idfilter, enabled, windowtitle)
            obj.figurehandle = fhandle;
            obj.slider = [];
            if nargin < 3
               idfilter = [];
            end
            if nargin < 4
               enabled = true; 
            end
            if nargin < 5
               windowtitle = ''; 
            end
           
            
            obj.idfilter = idfilter;
            obj.enabled = enabled;
            obj.windowtitle = windowtitle;
            
            obj.panel = uipanel(fhandle, 'DeleteFcn', @(o, e) obj.destructor(), 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            
            if isa(settings, 'Session')
                session = settings;
                obj.listener = session.addlistener('settingsChanged', @(src, ev) obj.settingsChanged(src));
                obj.locklistener = session.addlistener('lockChanged', @(src, ev) obj.lockChanged(session));
                obj.session = session;
                settings = session.settings;
            end
            
            
            
            obj.namemap = struct('cat0', '<no settings found>');
            categories = fieldnames(DefaultValues.SECTIONID);
            for k = 1:length(categories)
               obj.namemap.(sprintf('cat%i',DefaultValues.SECTIONID.(categories{k}))) = DefaultValues.SECTIONHEADER.(categories{k});
            end
            
            obj.setup(settings);
            
            scrollAmmount = DefaultValues.LETTERDIMENSION(2) * 1;
            if ishandle(fhandle)
            set(fhandle, 'Visible', 'on', 'WindowScrollWheelFcn', @(o, e) onScrollEvent(obj, e.VerticalScrollCount, scrollAmmount));
            end
           end
        
        function lockChanged(obj, session)
            if session.isLocked()
                set(allchild(obj.panel), 'Enable', 'off');
            else
                obj.settingsChanged(session);
            end
            
        end
        
        function destructor(obj)
           delete(obj.settinglistener);
           delete(obj.locklistener);
           delete(obj.listener);
           delete(obj); 
        end
        
        function settingsChanged(obj, session)
            delete(allchild(obj.panel));
            obj.setup(session.settings);
        end
        
        function setup(obj, settings)
            delete(obj.settinglistener); obj.settinglistener = [];
            catlist = settings.getSettings();
            subcatlist = [];
            catid = 0;
            
            for i = 1:length(catlist)
                if any(catlist{i}{1} == obj.idfilter) || isempty(obj.idfilter)
                    catid = catlist{i}{1};
                    titlename = catlist{i}{2}{1};
                    subcatlist = [subcatlist, catlist{i}(2:end)];
                end
            end
            grid = {};
            for i = 1:length(subcatlist)
                subcat = subcatlist{i}{1};
                
                mainid = settings.getSetting(subcatlist{i}{2}).getGroupID();
                dvlabel = sprintf('s%i_%s', mainid, subcat);
                if isfield(DefaultValues.SECTIONNAMES, dvlabel)
                   sectionname =  DefaultValues.SECTIONNAMES.(dvlabel);
                else
                    sectionname = dvlabel;
                end
                if ~isempty(sectionname)
                    grid{end+1,1} = GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10);
                    
                    sectionsettings = DefaultValues.SECTIONNAMESETTING;
                    grid{end+1,1} = uicontrol(obj.panel, 'style', 'text', 'String' , sectionname,   sectionsettings{:});
                end
                subgrid = {};
                
                for j = 2:length(subcatlist{i})
                    settingname = subcatlist{i}{j};                   
                    box = settings.renderGUI(obj.session, settingname, obj.panel, obj.options, obj.suggestions);
                    if ~isempty(box)
                        subgrid{end+1,1} = box;
                    end
                end
                
                sectionbox = LayoutBox(subgrid);
                %sectionbox = sectionbox.alignBoxes();
                grid{end+1,1} = sectionbox;
                
                
            end
            if ~isempty(obj.session)
                if catid == DefaultValues.SECTIONID.STARTER
                    warninghandle = uicontrol(obj.panel, 'style', 'text', 'String', '~~', 'ForegroundColor',[0.8,0, 0],'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR );
                    grid = [{warninghandle}; grid];
                    obj.settinglistener = settings.addlistener('settingChanged', @(o,e) GUIRenderSettings.onSettingChanged(obj.session, warninghandle));
                    GUIRenderSettings.onSettingChanged(obj.session, warninghandle)
                elseif catid == DefaultValues.SECTIONID.CONTINUER || catid == DefaultValues.SECTIONID.INTEGRATOR
                    branchmanager = obj.session.branchmanager;
                    [labels, names, name] = branchmanager.getConfVariants(obj.session.getCompConf().getSolutionLabel());
                    if ~isempty(labels)
                        
                        listbox = GUISelectCompConfListBox(obj.panel, obj.session, branchmanager, names, labels, 'HorizontalAlignment' , 'left', obj.options{:});
                        subgrid = {uicontrol( obj.panel, 'Style', 'text', 'String', name, 'HorizontalAlignment', 'left', obj.options{:},'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR ), listbox.handle };
                        box = LayoutBox(subgrid);
                        box.widths(1) = max(obj.suggestions.labelsize, box.widths(1));
                        grid = [{box}; grid];
                    end
                end
                
                
            
            end
            mainbox = LayoutBox(grid);
            obj.panel.UserData = mainbox;
            obj.panel.ResizeFcn = @(o, e) obj.onResize(o);
            
            if ~obj.enabled
               set(allchild(obj.panel), 'Enable', 'Inactive'); 
            end
            if isempty(obj.windowtitle)
                set(obj.figurehandle, 'Name', obj.namemap.(sprintf('cat%i', catid)));
            else
                set(obj.figurehandle, 'Name', obj.windowtitle);
            end
            
            if length(grid) == 0
                close(obj.figurehandle)
            else
                obj.onResize(obj.panel);
            end
        end
        
        function onResize(obj, source)
            box = source.UserData;
            delete(obj.slider);
            source.Units = 'pixels';
            obj.slider = box.doSliderLayout(source);
            source.Units = 'normalized';
            
            global slid; slid = obj.slider; %FIXME TODO REMOVE
        end
        
    end
    methods(Static)

        
        function onSettingChanged(session, handle)
            msg = session.computation.check(session.settings);
            if isempty(msg); msg = ' '; end
            
            set(handle, 'String', msg, 'TooltipString', msg);
            
        end
    end

end

function onScrollEvent(panel, direction, amount)
    %direction: 1 if scroll down, -1 if scroll up.
    if isempty(panel.slider); return; end
    
    newvalue = panel.slider.Value - direction*amount;
    newvalue = min(max(panel.slider.Min, newvalue), panel.slider.Max);
    panel.slider.Value = newvalue;
    panel.slider.Callback(panel.slider, []);
    

   
end

