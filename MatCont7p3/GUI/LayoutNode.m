classdef LayoutNode < handle
    
    properties
        wh
        ww
        index
        handles
        enabled = 1
        backup
        
        horizontal
        
        options = struct();
        
        slider = []
        
    end
    
    methods
        
        function destructor(obj)
            nrChildren = length(obj.handles);
            if ((nrChildren > 0) && ~obj.isHandleNode())
                children = vertcat(obj.handles{:});
                for i = 1:length(children)
                    children(i).destructor();
                end
            end
            delete(obj);
        end
        
        
        function obj = LayoutNode( wh , ww , varargin)
            obj.wh = wh;
            obj.ww = ww;
            obj.index = 0;
            obj.handles = cell(1,0);
            
            if (isempty(varargin))
                obj.horizontal = 1;
            else
                obj.horizontal = 0;
            end

            obj.options.minsize = [40,20];
            obj.options.halign = 'c';
            obj.options.sliderthreshold = [-Inf, -Inf];
            obj.options.margin = [0,0];
        end
        function setWeights(obj, wh , ww)
            obj.wh = wh;
            obj.ww = ww;
        end
        
        function addHandle(obj, wh , ww, handle , varargin)
            newnode = LayoutNode (wh, ww);
            newnode.index = 1;
            newnode.handles( newnode.index )  = {handle};
            newnode.setOptions(varargin{:});
            obj.addNode( newnode )
        end
        
        function setHandle(obj, handle, varargin)
            obj.index = 1;
            obj.handles( obj.index ) = {handle};
            obj.setOptions(varargin{:});
        end
        
        
        
        function addGUIobject(obj, wh, ww , guiobj , varargin)
           obj.addHandle( wh, ww , guiobj.handle , varargin{:}); 
        end
        
        function addNode(obj  , newnode)
           obj.index = obj.index + 1;
           obj.handles( obj.index ) = {newnode};
        end
        
        
        function makeLayoutHappen(obj , panelpos, callerid)
            delete(obj.slider); obj.slider = [];
           
           if ( panelpos(4) < obj.options.sliderthreshold(2))
              height = panelpos(4);
              panelpos(2) = height - obj.options.sliderthreshold(2);
              panelpos(4) = obj.options.sliderthreshold(2);
              panelpos(3) = panelpos(3) - 20;
              
              obj.slider(1) = LayoutSlider(obj.options.panel, [panelpos(3) , 0 , 20 , height] , obj, panelpos , (obj.options.sliderthreshold(2) - height) , 2 , (obj.options.sliderthreshold(2) - height)); 
           
           else
               panelpos(2) = 0;
           end
           
            if ( panelpos(3) < obj.options.sliderthreshold(1))
              width = panelpos(3);
              panelpos(1) = (width - obj.options.sliderthreshold(1));
              panelpos(3) = obj.options.sliderthreshold(1);
              panelpos(2) = panelpos(2) + 20;
              panelpos(4) = panelpos(4) - 20;
              obj.slider(end+1) = LayoutSlider(obj.options.panel, [0 , 0 , width, 20] , obj, panelpos , (obj.options.sliderthreshold(1) - width) , 1 , 0); 
              panelpos(1) = 0;
            else
              panelpos(1) = 0;
            end
            
            
            obj.distributePos(panelpos);
        end
        
        function distributePos(obj,pos)
           if (obj.horizontal)
                obj.distributePos_Horizontal(pos);
           else
                obj.distributePos_Vertical(pos);
           end
            
        end
            
        
        function distributePos_Vertical(obj , pos)
            if (obj.enabled == 0)
                return
            end
            
            if (isfield(obj.options, 'Add'))
                distributePos_Vertical_A(obj, pos)
                return;
            end
            
            left = pos(1); bottom = pos(2); width = pos(3); height = pos(4); 
            
            nrChildren = length(obj.handles);
            if (nrChildren == 0)
                return      
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                obj.setPosition(h,pos);
                
            else
            
                children = vertcat(obj.handles{:});
                whsum = sum(vertcat(children(:).wh));
                if (whsum == 0), return; end

                cLeft = left;
                cBottom = bottom;

                for i=length(children):-1:1
                   cWidth = width;
                   cHeight = height * ( children(i).wh / whsum);
                   children(i).distributePos( [cLeft cBottom cWidth cHeight]);
                   cBottom = cBottom + cHeight; 
                end
            end
            
        end
                
        function distributePos_Vertical_A(obj, pos)
 
            
         if (obj.enabled == 0)
                return
         end    
         left = pos(1); bottom = pos(2); width = pos(3); height = pos(4); 
         nrChildren = length(obj.handles);
            if (nrChildren == 0)
                return     
                
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                obj.setPosition(h,pos);
                
            else
            
                children = vertcat(obj.handles{:});
                

                cLeft = left ;
                cBottom = bottom + height;
                
                
                for i=1:length(children)
                   whsum = sum(vertcat(children(i:end).wh));
                   if (whsum == 0), return; end
                   cWidth = width;
                   cHeight = (cBottom - bottom) * ( children(i).wh / whsum);
                   [w,h] = getPrefSize(children(i));
                   
                   if (h < cHeight)
                       cHeight = h;
                   end
                   
                   cBottom = cBottom - cHeight; 
                  
                   children(i).distributePos( [cLeft cBottom cWidth cHeight]);
                   
                end
            end                
                
                
                
            
        end
        
        function distributePos_Horizontal(obj , pos)
            if (obj.enabled == 0)
                return
            end
            
            left = pos(1); bottom = pos(2); width = pos(3); height = pos(4);   
            
            nrChildren = length(obj.handles);
            
            if (nrChildren == 0)
                return
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                obj.setPosition(h,pos);
                
            else
                 children = vertcat(obj.handles{:});
                 wwsum = sum(vertcat(children(:).ww));
                 
                 cBottom = bottom;
                 cLeft = left;
                 
                 for i = 1:length(children)
                    cHeight = height;
                    cWidth = width * ( children(i).ww / wwsum);
                    children(i).distributePos(  [cLeft cBottom cWidth cHeight] );
                    cLeft = cLeft + cWidth;
                 end
                 
                
            end
        end
        

        
        
        function setEnabled(obj , bool)
            if (bool == obj.enabled)
                return;
            end
            obj.enabled = bool;
            
           if (bool) 
                obj.ww = obj.backup.ww;
                obj.wh = obj.backup.wh;
           else
               obj.backup.ww = obj.ww;
               obj.backup.wh = obj.wh;
               obj.ww = 0;
               obj.wh = 0; %REDO , hackish
           end
           
           
           nrChildren = length(obj.handles);
           if (nrChildren == 0)
               return
           elseif ((nrChildren == 1) && obj.isHandleNode())
               h = obj.handles{1}; %handle case
               set(h,'Visible' , bool2str( bool ));
               
           else
                children = vertcat(obj.handles{:});
                for i = 1:length(children)
                   children(i).setEnabled(bool); 
                end               
               
           end
           
        end
        
        function applyToHandles(obj, handlefunc)
            nrChildren = length(obj.handles);
            if (nrChildren == 0)
                return
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                type = get(h, 'Type');
                if (strcmp(type,'uipanel'))
                  data = get(h , 'UserData');
                  data.node.applyToHandles(handlefunc);
                else
                    handlefunc(h);
                end
                
            else
                children = vertcat(obj.handles{:});
                for i = 1:length(children)
                    children(i).applyToHandles(handlefunc);
                end
                
            end 
        end
        
        function makeEmpty(obj)

            for i = 1:length(obj.handles)
               applyToHandles(obj.handles{i}, @(h) delete(h)); 
               destructor(obj.handles{i});
            end
            obj.index = 0;
            obj.handles = {};
            
        end
        
        function bool = isHandleNode(obj)
            %{
            bool = ~isobject( obj.handles{1} );
            %}
            bool = ~isa(obj.handles{1}, 'LayoutNode');
            
        end
        
        function setPosition(obj,handle, pos)
            assert( obj.enabled ~= 0 );
            [w,h] = obj.getExtent(handle);
            
            
            width = min(pos(3) , w);
            height = min(pos(4) , h);
            
            bottom = pos(2) + (pos(4) - height)/2 ;
            if (obj.options.halign == 'l')
                left = pos(1);
            elseif(obj.options.halign == 'r')
                left = pos(1) + (pos(3) - width - 5);
            else
                left = pos(1) + (pos(3) - width)/2;
            end
            
            %add margins:
            left = max(pos(1) + obj.options.margin(1), left);
            width = min(pos(3) - 2*obj.options.margin(1) , width);
            bottom = max(pos(2) + obj.options.margin(2), bottom);
            height = min(pos(4) - 2*obj.options.margin(2) , height);
            
            set(handle ,'Position' ,  [left,bottom,width,height]);
        end
        
        function [w,h] = getExtent(obj,handle)
            try
                type = get(handle, 'Type');
                if (strcmp(type,'uicontrol') || strcmp(type,'uitable'))
                
                    ext = get(handle, 'Extent');
                    w = ext(3) * 1.05;
                    h = ext(4) * 1.05;
                    w = max(w , obj.options.minsize(1));
                    h = max(h , obj.options.minsize(2));
                elseif(strcmp(type,'uipanel'))
                    data = get(handle, 'UserData');
                    [w,h] = data.extentfunction();
                    w = w * 1.05;
                    h = h * 1.05;
                    
                    w = max(w , obj.options.minsize(1));
                    h = max(h , obj.options.minsize(2));                    
                    
                else
                    %disp('Unknown UI object encountered');
                    w=Inf; h = Inf;
                end
            catch
                %disp('unExtentable object encountered');
                w=Inf; h = Inf;
            end
        end

         function [w,h] = getTrueExtent(obj,handle)
            try
                type = get(handle, 'Type');
                if (strcmp(type,'uicontrol') || strcmp(type,'uitable'))
                
                    ext = get(handle, 'Extent');
                    w = ext(3) * 1.05;
                    h = ext(4) * 1.05;
                elseif(strcmp(type,'uipanel'))
                    data = get(handle, 'UserData');
                    [w,h] = data.extentfunction();
                    w = w * 1.05;
                    h = h * 1.05;
                    
                else
                    %disp('Unknown UI object encountered');
                    w=0; h = 0;
                end
            catch err
                %disp('unExtentable object encountered');
                w=0; h = 0;
            end
        end       
        
        
        
        function setOptions(obj,varargin) %good citizen
            for i = 1:2:length(varargin)
                obj.options.(varargin{i}) = varargin{i+1};
            end
        end
        function disableFromBelow(obj)
           enabled = 0;
           for i = 1:length(obj.handles)
               if (obj.handles{i}.enabled == 1)
                   enabled = 1;
               end
           end
           obj.setEnabled(enabled);
        end
        
        function [w,h] = getPrefSize(obj) %experimental, flawed
            nrChildren = length(obj.handles);
            if ((nrChildren == 0) || (obj.enabled == 0))
                w=0; h=0;
                return
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                [w,h] = obj.getExtent(h);

            else
                children = vertcat(obj.handles{:});
                w=0;h=0;
                for i = 1:length(children)
                    [w1,h1] = children(i).getPrefSize();
                    
                    if (obj.horizontal)
                        w = w + w1;
                        h = max(h , h1);
                    else
                        h = h + h1;
                        w = max(w , w1);
                    end
                end
                
            end 
            
        end
        function [w,h] = getTruePrefSize(obj)
            
            nrChildren = length(obj.handles);
            if ((nrChildren == 0) || (obj.enabled == 0))
                w=0; h=0;
                return
            elseif ((nrChildren == 1) && obj.isHandleNode())
                h = obj.handles{1}; %handle case
                [w,h] = obj.getTrueExtent(h);

            else
                children = vertcat(obj.handles{:});
                w=0;h=0;
                for i = 1:length(children)
                    [w1,h1] = children(i).getTruePrefSize();
                    
                    if (obj.horizontal)
                        w = w + w1;
                        h = max(h , h1);
                    else
                        h = h + h1;
                        w = max(w , w1);
                    end
                end
                
            end 
 
        end
        
    end
    methods(Static)
        function normalize(handle)
            set( findall(handle) , 'Units' , 'normalized');

        end
    end

end
         
function result = bool2str( bool)
if  (bool)
    result = 'on';
else
    result = 'off';
end
end

function h = LayoutSlider(parent , sliderpos , mainnode , panelposinit , maxval , index , initval)
min_step = min(0.5 , max(20/maxval,0.01));
max_step = min(0.7 , max(10*min_step , 0.1));
mainnode.options.pos = panelposinit;

 h = uicontrol(parent , 'Style' , 'slider' , 'Position', sliderpos , 'Max' , maxval , 'Min' , 0 , 'Value' , initval ...
    ,'SliderStep' , [min_step max_step],  'Callback' , @(ho,e)  sliderCallBack(ho, mainnode , panelposinit , index, maxval) );
end

function sliderCallBack(sliderhandle  , mainnode , panelposinit,  index , maxval)

val = get(sliderhandle , 'Value');

mainnode.options.pos(index) = panelposinit(index) + (maxval - val);
mainnode.distributePos(mainnode.options.pos);

end

