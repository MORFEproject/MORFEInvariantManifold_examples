classdef GUIPanel < handle

    properties
        handle
        mainnode
        borderwidth
        borderheight
        bannerwidth
    end
    
    methods
        function obj = GUIPanel(parent,bannerw , border, varargin)
            obj.handle = uipanel(parent, 'Unit' , 'Pixel' , 'DeleteFcn' , @(o,e) obj.destructor() , 'ResizeFcn' , @(o,e) obj.doLayout() , varargin{:});
            obj.mainnode = LayoutNode(-1,-1,'vertical'); 
            obj.bannerwidth = bannerw;
            
            if (length(border) == 1)
               obj.borderwidth = border;
               obj.borderheight = border;
            else
               obj.borderwidth = border(1);
               obj.borderheight = border(2);
            end
            
            data.extentfunction = @() obj.getExtent();
            data.node = obj.mainnode;
            
            set(obj.handle , 'UserData' , data);
        end
        
        
        
        function doLayout(obj)
            pos = get(obj.handle, 'Position');
            bw = max( get(obj.handle, 'BorderWidth') , obj.borderwidth) ;
            bh = max( get(obj.handle, 'BorderWidth') , obj.borderheight) ;
            if (~isempty(pos))
                pos(4) = max(pos(4) - obj.bannerwidth,1);
                pos(1) = bw;
                pos(2) = bh;
                pos(3) = max(pos(3) - 2*bw,1);
                pos(4) = max(pos(4) - 2*bh,1);
                obj.mainnode.distributePos(pos);
            end
        end
        
        function destructor(obj)
           delete(obj.mainnode);
        end
        
        function [w,h] = getExtent(obj)
            [w,h] = obj.mainnode.getPrefSize();
            h = (h + obj.bannerwidth + 2* obj.borderheight);
            w = (w + 2*obj.borderwidth) ;
            
            
        end
    end
    
end

