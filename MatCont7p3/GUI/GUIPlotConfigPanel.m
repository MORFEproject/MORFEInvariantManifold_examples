classdef GUIPlotConfigPanel < handle

    
    properties
        dimension;
        panelhandle;
        plotconf;
        regionLeftHandle;
        regionRightHandle;
        

        mainnode
        
        catpopuphandle
        subcatpopuphandle
        modifierpopuphandle
        
        selectionListener;
    end
    
    
    methods
        function obj = GUIPlotConfigPanel(parent , plotconf , dim  ,varargin)
            obj.plotconf = plotconf;
            obj.dimension = dim;
            obj.panelhandle = uipanel(parent, 'Unit' , 'Pixel' ,'Title' , plotconf.getAxTitle(dim)  ,...
                'DeleteFcn' , @(o,e) obj.destructor() , varargin{:});
            
            mainnode = LayoutNode(-1,-1, 'vertical');
            obj.mainnode = mainnode;
            
            subnode = LayoutNode(1,0);
            
            obj.catpopuphandle = uicontrol(obj.panelhandle , 'Units', 'Pixels' , 'Style' , 'popupmenu' ...
                , 'Callback' , @(o,e) obj.catPopupCallback(o,e)   );
            obj.subcatpopuphandle = uicontrol(obj.panelhandle , 'Units', 'Pixels' , 'Style' , 'popupmenu' ...
                , 'Callback' , @(o,e) obj.sectionPopupCallback(o,e)   );
            obj.modifierpopuphandle = uicontrol(obj.panelhandle , 'Units', 'Pixels' , 'Style' , 'popupmenu' ...
                , 'Callback' , @(o,e) obj.sectionPopupCallback(o,e)   );
            

            subnode.addHandle(0,2, obj.subcatpopuphandle,   'minsize' , [Inf 20] , 'margin', [5 5]);	
            subnode.addHandle(0,2, obj.modifierpopuphandle, 'minsize' , [Inf 20] , 'margin', [5 5]);
            subnode.addHandle(0,4, obj.catpopuphandle,      'minsize' , [Inf 20] , 'margin', [5 5]);
                       
            mainnode.addNode(subnode);
            subnode = LayoutNode(1,0);
            subnode.addHandle(0,1,uicontrol(obj.panelhandle, 'Units','Pixels', 'Style' , 'text' , 'String' , 'Range:'));
            
            obj.regionLeftHandle = uicontrol(obj.panelhandle, 'Units','Pixels', 'Style'  , 'edit' , 'String' , 'a'  , 'Callback' ,@(o,e) obj.regionChangeCallback(plotconf) );
            obj.regionRightHandle = uicontrol(obj.panelhandle, 'Units','Pixels', 'Style'  , 'edit' , 'String' , 'b'   , 'Callback' ,@(o,e) obj.regionChangeCallback(plotconf));
            
            subnode.addHandle(0,2,obj.regionLeftHandle, 'minsize' , [Inf 20] , 'margin' , [5,5]);
            subnode.addHandle(0,1,uicontrol(obj.panelhandle,'Units','Pixels', 'Style' , 'text' , 'String' , '...' , 'Fontsize' ,13));
            subnode.addHandle(0,2,obj.regionRightHandle, 'minsize' , [Inf 20], 'margin' , [5,5]);
            
            mainnode.addNode(subnode);
            
            data.extentfunction = @() obj.getExtent();
            data.node = obj.mainnode;
            set(obj.panelhandle , 'UserData' , data);
            
            
            set(obj.panelhandle, 'ResizeFcn' , @(o,e) obj.doLayout());
            obj.doLayout();
            
            obj.configurePopups();
            obj.selectionListener = obj.plotconf.addlistener('selectionChanged', @(o, e) obj.configurePopups());
        end
        
        function configurePopups(obj)
            [l,r] = obj.plotconf.getRegion(obj.dimension);
            set(obj.regionLeftHandle, 'String', num2str(l, '%.8g'));
            set(obj.regionRightHandle, 'String', num2str(r, '%.8g'));
            
            set([obj.subcatpopuphandle, obj.modifierpopuphandle], 'Visible', 'off');
            [catlist, selection] = obj.plotconf.getCategories(obj.dimension);
            
            if isempty(selection) 
               catlist = [{'', ''}; catlist];
               catindex = 1;
            else
               catindex = find(strcmp(catlist(:, 1), selection));

               [subcatlist, selection, singleton] = obj.plotconf.getSubcatItems(obj.dimension);
               subcatindex = find(strcmp(subcatlist(:, 1), selection));
               set(obj.subcatpopuphandle, 'Value', subcatindex, 'String', subcatlist(:, 2), 'visible', 'on', 'UserData', subcatlist(:,1) , 'visible', CLbool2text(~singleton));
               
               [itemlist, selection] = obj.plotconf.getItems(obj.dimension);
               newitemlist = itemlist;
               for k = 1:length(itemlist)
                  newitemlist{k} = replace(newitemlist{k}, '_', ' '); 
               end
               itemindex = find(strcmp(itemlist, selection));
               set(obj.modifierpopuphandle, 'Value', itemindex, 'String', newitemlist, 'UserData', itemlist, 'visible', CLbool2text(length(itemlist) > 1 ));
            end
            set(obj.catpopuphandle, 'Value', catindex, 'String', catlist(:, 2), 'UserData', catlist(:, 1));
            
        end
        
        function catPopupCallback(obj, ~, ~)
            index = get(obj.catpopuphandle, 'value');
            catname = obj.catpopuphandle.UserData{index};
            obj.plotconf.setCategory(obj.dimension, catname);
           
        end
        
        
        function sectionPopupCallback(obj,handle,e)
            index = get(obj.subcatpopuphandle, 'value');
            subcatname = obj.subcatpopuphandle.UserData{index};
            
            index = get(obj.modifierpopuphandle, 'value');
            item = obj.modifierpopuphandle.UserData{index};
            
            obj.plotconf.setSelection(obj.dimension, subcatname, item);
            cla(obj.plotconf.axeshandle);  % clear plot
        end
        
        
        function destructor(obj)
           %fprintf(2, 'destructor GPCP\n');
           delete(obj.mainnode); 
           delete(obj.selectionListener);
        end


        
        function [w,h] = getExtent(obj)
            [w,h] = obj.mainnode.getPrefSize();
            h = h + 20;
            w = w + 4;
        end
        
        
        function doLayout(obj)
            units = get(obj.panelhandle, 'Units');
          
            set(obj.panelhandle,'Units', 'Pixels');
            pos = get(obj.panelhandle, 'Position');
            if ~isempty(pos)
                pos(4) = max(pos(4) - 10,1);
                obj.mainnode.makeLayoutHappen(pos);
            end
            set(obj.panelhandle,'Units',units);
        end
        
        function regionChangeCallback(obj, plotconf)
            plotconf.setAutoFit(false);
            try
                region = [str2double(get(obj.regionLeftHandle, 'String')), str2double(get(obj.regionRightHandle, 'String'))];
                if diff(region) > 0
                    plotconf.setRegion(obj.dimension, [str2double(get(obj.regionLeftHandle, 'String')), str2double(get(obj.regionRightHandle, 'String'))]);
                end
            catch
            end
            
        end
    end

    
end

                

