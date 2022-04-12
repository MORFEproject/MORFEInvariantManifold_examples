classdef GUILineConfigPanel < handle

    properties
        handle;
        lineoptions;
        slider = [];
        mainbox;
    end
    
    events
       settingChanged
    end
    

    methods
        function obj = GUILineConfigPanel(parent , session , varargin)
            obj.handle = uipanel(parent, 'Unit' , 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, varargin{:});
            uimenu(parent, 'Label' , 'Help' , 'Callback' , @(o,e) obj.displayHelp());
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor() , 'ResizeFcn' , @(o,e) obj.onResize(o));
            obj.slider = [];
            
           
            obj.lineoptions = session.globalsettings.lineoptions;
            
            labels = obj.lineoptions.getLabels();
            
            grid = cell(2*length(labels), 1);
            
            separator = '';
            for k = 1:length(labels)
               label = labels{k};
               [~, name] = obj.lineoptions.getSetting(label);
               
               grid{2*k - 1} = uicontrol( obj.handle  , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , 'HorizontalAlignment' , 'left' , ...
                    'String' , [separator , name , ':'] );
               
                subgrid = cell(1, 2);
                subgrid{1} = uicontrol(obj.handle , 'Style' , 'pushbutton' , 'String' , 'Default','Callback' , @(o,e) obj.setDefault(label), 'TooltipString', ['restore the default value for "' name '"']);
                subgrid{2} = GUIEvalPlotOpsBox( obj.handle , ...
                              @(x) obj.lineoptions.setSetting(label , x) ...
                             , @() obj.lineoptions.getSetting(label) ...
                             , obj, 'settingChanged' );
                grid{2*k} = LayoutBox(subgrid);
                separator = newline; %\n
            end
            obj.mainbox = LayoutBox(grid);
            obj.mainbox.sourcepanel = obj;
            
            obj.mainbox.makeLayoutHappen( get(obj.handle , 'Position'))
            
            set(obj.handle, 'Units', 'normalized');
            scrollAmmount = DefaultValues.LETTERDIMENSION(2) * 1;
            set(parent, 'WindowScrollWheelFcn', @(o, e) onScrollEvent(obj, e.VerticalScrollCount, scrollAmmount));
            
        end
        
        
        function onResize(obj, handle)
           set(handle,'Units' , 'Pixels');
           pos = get(handle, 'Position');
           if (~isempty(pos)) 
                obj.mainbox.makeLayoutHappen( get(handle, 'Position'));
           end
           set(handle,'Units' , 'normalize');
        end
           
        
        
        function destructor(obj)
           delete(obj);
        end
        
        function setDefault(obj, label)
           obj.lineoptions.setDefault(label);
           obj.notify('settingChanged');
            
        end
        
    end
    methods(Static)
        function displayHelp()
           helpdlg( ...
             [ 'The options in these editboxes will be passed along to the MATLAB plot instructions.', ... 
             'The syntax of these instructions are the same as the one MATLAB uses for the instrunctions "plot", "line" and "text" and can be consulted in the documentation pages of these instructions.' ]...
           , 'Plot Settings');
            
        end
        
        
    end

end

function onScrollEvent(panel, direction, amount)
    %direction: 1 if scroll down, -1 if scroll up.
    if isempty(panel.slider) || ~isvalid(panel.slider); return; end
    
    slider = panel.slider.handle;
    newvalue = slider.Value - direction*amount;
    newvalue = min(max(slider.Min, newvalue), slider.Max);
    slider.Value = newvalue;
    slider.Callback(slider, []);
    

   
end
