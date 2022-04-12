classdef LayoutBox
    properties
        grid;
        widths;
        heights;
        
        fillLastColumnHorizontal = true; 
        fillLastAlignRight = false;
        
        margins = [0, 0];
        margins2 = [0, 0];
        
        sourcepanel;
    end
    
    methods
        function obj = LayoutBox(grid)
            
            obj.grid = grid;
            [~, w, h] = obj.Extent();
            obj.widths = w + obj.margins(1);
            obj.heights = h + obj.margins(2);
            

        end
 
        function [ext, width, height] = Extent(obj)
            
            [n, m] = size(obj.grid);
            width = zeros(1, size(obj.grid,2));
            height = zeros(size(obj.grid,1), 1);
            for i = 1:n
                for j = 1:m
                    if ~isempty(obj.grid{i,j})
                        ext = obj.grid{i,j}.Extent;
                        width(j) = max(width(j), ext(3));
                        height(i) = max(height(i), ext(4));
                    end
                end
            end
            
            ext = [0, 0, sum(width), sum(height)];
        end
        
        
        function setPosition(obj, position)
            x = position(1); y = position(2);
            maxwidth = position(3);
            

            [n, m] = size(obj.grid);
            
            for i = 1:n
                xt = x;
                y = y - obj.heights(i) - obj.margins2(2);
                for j = 1:m
                    if ~isempty(obj.grid{i,j})

                        if j == m && obj.fillLastColumnHorizontal
                            
                           if ~obj.fillLastAlignRight
                                w = min(max(obj.widths(j), maxwidth - xt - x), maxwidth);
                           else
                                w = obj.widths(j); 
                                xt = max(xt, maxwidth - (w + 1) - DefaultValues.LETTERDIMENSION(2)*2);
                           end
                        else
                           w = obj.widths(j); 
                        end
                        
                        set(obj.grid{i,j}, 'Position', [xt, y, w, obj.heights(i)]);
                    end
                    xt = xt + obj.widths(j);
                end
                
            end
            
        end
        
        function sliderhandle = doSliderLayout(obj,panelhandle)
            panelpos = get(panelhandle, 'Position');
            sliderhandle = [];
            preferredheight = sum(obj.heights);
            if (panelpos(4) < preferredheight)
                slider = GUISlider(obj, panelhandle, panelpos, preferredheight);
                sliderhandle = slider.handle;
            else
                set(obj, 'Position', panelpos);
            end
        end
        
        function set(obj, varargin)
            ii = find(strcmp(varargin, 'Position'));
            if ii
                position = varargin{ii+1};
                position(2) = position(2) + position(4);
                obj.setPosition(position);
            end
        end
        function makeLayoutHappen(obj , panelpos)
            panelpos(4) = floor(panelpos(4)*0.99);
            
            if ~isempty(obj.sourcepanel.slider) && isvalid(obj.sourcepanel.slider)
                delete(obj.sourcepanel.slider.handle);
            end
            preferredheight = sum(obj.heights);
            panelpos = [1 1 panelpos(3:4)];
            
            if (panelpos(4) < preferredheight)
                obj.sourcepanel.slider = GUISlider(obj,obj.sourcepanel.handle , panelpos, preferredheight);
                
            else
                obj.sourcepanel.slider = [];
                set(obj, 'Position', panelpos);
                
            end
            
        end
        function destructor(obj)
            
        end

    end
end

