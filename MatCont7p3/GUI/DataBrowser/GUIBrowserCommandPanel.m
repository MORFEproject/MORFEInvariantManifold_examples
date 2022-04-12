classdef GUIBrowserCommandPanel
    
    properties
        handle
        model        
        eventlistener;
    end
    
    methods
        
        function obj = GUIBrowserCommandPanel(parent, model, varargin)
            obj.model = model;
            obj.handle = uipanel(parent, 'Unit' , 'Pixels'  , varargin{:},'ResizeFcn' ,  @(o,e) obj.reconfigure());    
            obj.eventlistener = model.addlistener('listChanged' , @(o,e) obj.reconfigure());
            
            obj.configure(get(obj.handle , 'Position' ));
            
        end
        
        function reconfigure(obj)
           if ~isempty(obj.handle); delete(allchild(obj.handle)); end

           units = get(obj.handle, 'Unit');
           set(obj.handle, 'Unit' , 'pixel');
           pos = get(obj.handle , 'Position');
           if (~isempty(pos))
                obj.configure(pos);
           end
           set(obj.handle, 'Unit' , units);
        end
        
        function configure(obj, pos)
            ops_ = obj.model.current.getOperations();
            ops_inh = obj.model.current.getInheritedOperations();
            ops = [ops_ ops_inh ];
            
             if (obj.model.getNrElements() > 0)
                ops_item = obj.model.current.getItemOperations();
             else
                 ops_item = [];
             end
          
             ops_len = length(ops) + length(ops_item);
             max_width = pos(3) / ops_len;
             panelheight = min(pos(4) , 23);
             
            left = 0;
            
            for i = 1:length(ops)
                s = ops{i};
                pb = GUIPushButton(obj.handle, s.cmd , s.label , obj.model , @() 1 ...
                    , 'indexChanged');
                
               ext = get(pb.handle , 'Extent');
               width = max( min(max_width , (ext(3)+5) ) , 80);
               set(pb.handle , 'Position' , [left 0 width panelheight]);
               
               left = left + width;
                
            end
            %FIXME TODO:  merge code
            
            for i = 1:length(ops_item)
                s = ops_item{i};
                
                
                pb = GUIPushButton(obj.handle, s.cmd , s.label , obj.model , @() obj.model.current.isValidSelection() ...
                    , 'indexChanged');
                ext = get(pb.handle , 'Extent');
                width = max( min(max_width , (ext(3)+5) ) , 80);
                set(pb.handle , 'Position' , [left 0 width panelheight]);
                left = left + width;
            end
        end
        
             
    end
    
end

