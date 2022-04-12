classdef GUIBrowserList < handle
    
    properties
        handle
        prevtime
        previndex;
        
        browsermodel;
        
        eventlisteners = {};
    end
    
    methods
        function obj = GUIBrowserList(parent, browsermodel , varargin)
            
            obj.handle = uicontrol(parent, 'Style' , 'listbox' ,  'Unit' , 'Pixels' , ... 
            'Callback' , @(o,e) obj.selectCallback(),'DeleteFcn' , @(o,e) obj.destructor(),varargin{:},...
            'KeyPressFcn' , @(o,e) keypress(o,e,browsermodel));
      

            obj.setList(browsermodel.getSelectList(), browsermodel.getSelectedIndex());
            set(obj.handle , 'TooltipString' , browsermodel.getCurrentTooltipString() );
                
            obj.eventlisteners{end+1} = browsermodel.addlistener('listChanged', @(o,e) obj.listChanged());
            obj.eventlisteners{end+1} = browsermodel.addlistener('indexChanged' ,@(o,e) obj.indexChanged());
        
            obj.prevtime = clock;
            obj.previndex = -9;
            obj.browsermodel = browsermodel;
        end

        
        
        %%%%%%%
        function setList(obj,list, index) 
            displist = cell(length(list)+2, 1);
            displist{1} = '..';
            displist{2} = '  ';
            
           for i = 1:length(list)
                displist{i+2} = deblank(list{i}); %list could contain nasty linefeeds and newlines.
           end
           set(obj.handle, 'String' , displist);
           
           if (isempty(list))
                set(obj.handle, 'Value', 1);
                
           else
               obj.setIndex(index);
           end
     
           

        end
        function listm = getList(obj)
           list = get(obj.handle, 'String');
           listm = list(3:end);
        end
        function indexm = getIndex(obj)
            indexm = get(obj.handle,'Value') - 2;
        end
        function setIndex(obj,index)
            set(obj.handle, 'Value' , index + 2);
        end
        %%%%%%%%
        function selectCallback(obj)
           index = obj.getIndex();
           time = clock;
           
           
           if (index == obj.previndex) 
               if (etime(time,obj.prevtime) <= 0.5)
                    obj.previndex = -9;
                    if (index >= 1)
                        obj.browsermodel.selectItem(index);
                    elseif (index == -1)
                        obj.browsermodel.goUp();
                    end
                    return %! belangrijk
                end
           elseif (index >= 1)
              obj.browsermodel.setSelectedIndex(index);
           else
              %obj.browsermodel.unSetIndex(); 
           end
           
           obj.prevtime = time;
           obj.previndex = index;
        end
        
        
        function listChanged(obj)
           obj.setList(obj.browsermodel.getSelectList(), obj.browsermodel.getSelectedIndex());
           set(obj.handle , 'TooltipString' , obj.browsermodel.getCurrentTooltipString() );
        end
        
        function indexChanged(obj)
            index = obj.browsermodel.getSelectedIndex();
            index_list = obj.getIndex();
            if (index_list ~= index)
                obj.setIndex(obj.browsermodel.getSelectedIndex());
            end
        end
        
        function destructor(obj)
            for i = 1:length(obj.eventlisteners)
               delete(obj.eventlisteners{i}); 
            end
            
            delete(obj);
        end

                
    end
    
end
function keypress(handle, event, browsermodel)
browsermodel.keypress(event);


end
