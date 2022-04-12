classdef GUINumericPanel < handle
    properties

        panel;
        model;
        layouthandle;
        setters = {};
        
        %options =  {'Units', 'Pixels', 'Visible', 'on'};
        %suggestions = struct('labelsize', 150 );
        %widthmargin = 10;
        %session = [];
        
        
        eventlistener = [];
        eventlistener2 = [];
        
        slider;
        
    end
    
    
    
    methods
        function obj = GUINumericPanel(parent, numericmodel)
            obj.panel = uipanel(parent, 'DeleteFcn' , @(o,e) obj.destructor(), 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            obj.layouthandle = uimenu(parent, 'Label', 'Layout');
            obj.changeNumericModel(numericmodel);
            obj.slider = [];
            

            scrollAmmount = DefaultValues.LETTERDIMENSION(2) * 1;
            set(parent, 'WindowScrollWheelFcn', @(o, e) onScrollEvent(obj, e.VerticalScrollCount, scrollAmmount));

        end
        
        function changeNumericModel(obj, numericmodel)
            obj.model = numericmodel;
            delete(obj.eventlistener);
            delete(obj.eventlistener2);
            obj.eventlistener = obj.model.addlistener('visibilityChanged' , @(srv,ev) obj.setupPanel());
            obj.eventlistener2 = obj.model.addlistener('categoriesChanged' , @(srv,ev) obj.setup());
            obj.setup();
        end
        
        function destructor(obj)
           delete(obj.layouthandle);
           delete(obj.eventlistener);
           delete(obj.eventlistener2);
           delete(obj);
           
        end
        
        
        function setup(obj)
            obj.setupPanel();
            obj.setupLayoutMenu();
        end
        
        function setupLayoutMenu(obj)
            delete(allchild(obj.layouthandle));
            catnames = obj.model.getCategories();
            for i = 1:length(catnames)
                GUISwitchMenuItem(obj.layouthandle, obj.model.getFullName(catnames{i}), @(b) obj.model.setCatVisible(catnames{i}, b), @() obj.model.isCatVisible(catnames{i}), obj.model, 'visibilityChanged');
            end
            
        end
        
        function setupPanel(obj)
            delete(allchild(obj.panel));
          
            obj.setters = {};
            catnames = obj.model.getCategories();
            
            grid = cell(0, 2);
            for i = 1:length(catnames)
                if obj.model.isCatVisible(catnames{i})
                    labels = obj.model.getLabels(catnames{i});
                    grid{end+1,2} = GUIEmptySpace(round(DefaultValues.LETTERDIMENSION(1)/2), 10);
                    grid{end+1,2} = uicontrol(obj.panel, 'style', 'text', 'String' , obj.model.getFullName(catnames{i}), 'horizontalalignment' , 'center', 'FontAngle', 'italic', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                    
                    for j = 1:size(labels, 1)
                        datalabel = uicontrol(obj.panel, 'style', 'text', 'String' , '-', 'horizontalalignment' , 'center', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
                        obj.setters{end+1} = @(data,s,ind,i) set(datalabel, 'String', sprintf('%.15g', labels{j, 2}(data{:},s,ind,i)));
                        grid(end+1, :) = {uicontrol(obj.panel, 'style', 'text', 'String' , labels{j, 1}, 'horizontalalignment' , 'left', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR) , datalabel};
                    end
                end
                
                
            end
            mainbox = LayoutBox(grid);
            obj.panel.UserData = mainbox;
            obj.panel.ResizeFcn = @(o, e) obj.onResize(o);
            obj.onResize(obj.panel);
        end
        
        
        function output(obj, data, s, ind)
            for i = 1:length(obj.setters)
                obj.setters{i}(data, s, ind, ind(end));
            end
        end
        
        function outputPoint(~, varargin)  %not implemented for Numeric.
        end    
        
        function b = isValid(obj)
            b = 1; %check if anything selected?
        end
        function onResize(obj, source)
            box = source.UserData;
            source.Units = 'pixels';
            delete(obj.slider);
            obj.slider = box.doSliderLayout(source);
            source.Units = 'normalized';
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
