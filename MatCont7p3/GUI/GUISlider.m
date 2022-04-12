classdef GUISlider < handle
    
    
    properties
        position
        layoutbox
        handle
    end
    
    methods
        function obj = GUISlider(layoutbox, parent, panelpos, preferredheight)
            
            
            height = panelpos(4);
            panelpos(2) = height - preferredheight;
            panelpos(4) = preferredheight;
            panelpos(3) = panelpos(3) - 20;
            
            
            sliderpos = [panelpos(3) , 0 , 20 , height];
            panelposinit = panelpos;
            maxval = (preferredheight - height);
            index = 2;
            initval = (preferredheight - height);
            
            min_step = min(0.5 , max(20/maxval,0.01));
            max_step = min(0.7 , max(10*min_step , 0.1));
            
            obj.handle = uicontrol(parent , 'Style' , 'slider' , 'Position', sliderpos , 'Max' , maxval , 'Min' , 0 , 'Value' , initval ...
                ,'SliderStep' , [min_step max_step],  'Callback' , @(h,e)  obj.sliderCallBack(h , panelposinit , index, maxval)...
                ,'DeleteFcn', @(o, e) delete(obj));
            
            
            obj.position = panelposinit;
            obj.layoutbox = layoutbox;
            
            
            set(layoutbox, 'Position', panelpos);
        end
        
        
        
        
        function sliderCallBack(obj, sliderhandle  , panelposinit,  index , maxval)
            val = get(sliderhandle , 'Value');
            obj.position(index) = panelposinit(index) + (maxval - val);
            set(obj.layoutbox, 'Position', obj.position);
            
        end
        
        
    end
    
end