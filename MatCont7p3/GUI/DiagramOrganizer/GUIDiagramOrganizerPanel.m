classdef GUIDiagramOrganizerPanel

    
    properties
        handle
        commodel
        
        %%%
        
        diagramlist1
        diagramlist2
        
        curvelist1
        curvelist2
    end
    
    methods
        function obj = GUIDiagramOrganizerPanel(parent , commodel , varargin)
            obj.commodel = commodel;
            obj.handle = uipanel(parent , 'Unit' , 'Pixels' , 'DeleteFcn' , @(o,e) obj.destructor() , varargin{:});
            
            
            
            obj.diagramlist1 = uicontrol(obj.handle , 'Style' , 'popup' , 'Unit' , 'Pixels', 'String' , {' '} , 'Value' , 1 , 'UserData' , 1 ... 
                ,'Callback' ,@(o,e) diagramSelectCallback(o,commodel)  );
            
            obj.diagramlist2 = uicontrol(obj.handle , 'Style' , 'popup' , 'Unit' , 'Pixels' ,'String' ,  {' '}  , 'Value' , 1 , 'UserData' , 2 ...
                ,'Callback' ,@(o,e) diagramSelectCallback(o,commodel)  );
           
            %%%%%%%% 
            
            obj.curvelist1 = uicontrol(obj.handle , 'Style' , 'listbox' , 'Unit' , 'Pixels', 'String' , {' '} , 'Value' , 1 , 'UserData' , 1 ...
                ,'Callback' ,@(o,e) curveSelectCallback(o,commodel) , 'KeyPressFcn' , @(o,e)keypress(o,e,commodel , 1) );
            obj.curvelist2 = uicontrol(obj.handle , 'Style' , 'listbox' , 'Unit' , 'Pixels' , 'String' , {' '} , 'Value' , 1 , 'UserData' , 2 ...
                ,'Callback' ,@(o,e) curveSelectCallback(o,commodel) , 'KeyPressFcn' , @(o,e)keypress(o,e,commodel , 2) );
            
            pb1 = GUIPushButton(obj.handle,  @(o,e) commodel.moveFromLeft()  , '>>' , commodel , @()  commodel.moveFromLeftValid()  , 'diagramSelectionChanged'  );
            pb2 = GUIPushButton(obj.handle,  @(o,e) commodel.moveFromRight() , '<<' , commodel , @()  commodel.moveFromRightValid() , 'diagramSelectionChanged'  );
           
            %%%%%%%% 
                        
            obj.configure(commodel);
            
            commodel.addlistener('diagramSelectionChanged' , @(ob,ev) diagramSelectionChanged(obj.curvelist1 , ob , ev));
            commodel.addlistener('diagramSelectionChanged' , @(ob,ev) diagramSelectionChanged(obj.curvelist2 , ob , ev));
            
            setListContextMenu(obj.curvelist1 , commodel , 1);
            setListContextMenu(obj.curvelist2 , commodel , 2);
            setDiagramContextMenu(obj.diagramlist1 , commodel , 1);
            setDiagramContextMenu(obj.diagramlist2 , commodel , 2);
            
            commodel.addlistener('diagramListChanged' , @(o,e) obj.configure(commodel));   
            
            mainnode = LayoutNode(1,1);
            left = LayoutNode(1,4,'vertical');
            middle = LayoutNode(1,1,'Vertical');
            right = LayoutNode(1,4,'Vertical');
            
            left.addHandle(1,0, obj.diagramlist1 , 'minsize' , [Inf,Inf]);
            left.addHandle(9,0 , obj.curvelist1 , 'minsize' , [Inf,Inf]);
            
            middle.addNode(LayoutNode(1,0))
            middle_ = LayoutNode(9,0,'Vertical');
            middle_.addGUIobject(1,1,pb1);
            middle_.addGUIobject(1,1,pb2);
            middle.addNode(middle_);
            
            middle.addGUIobject(1,0,pb1);
            middle.addGUIobject(1,0,pb2);
            
            right.addHandle(1,0,obj.diagramlist2 , 'minsize' , [Inf,Inf]);
            right.addHandle(9,0,obj.curvelist2 , 'minsize' , [Inf,Inf]);
            mainnode.addNode(left);
            mainnode.addNode(middle);
            mainnode.addNode(right);
            
            mainnode.makeLayoutHappen(get(obj.handle, 'Position'));
            mainnode.destructor();

            set( obj.handle , 'Units','normalize');
            set( allchild(obj.handle) , 'Units','normalize');
        end
        
        
        function configure(obj,commodel)
                syncDiagramList(obj, commodel);    
                syncCurveList(obj , commodel);
        end
                
        
        function destructor(obj)
        end
          
        function syncDiagramList(obj , commodel)
            
            diagrams = commodel.getDiagramList();
            set([obj.diagramlist1 obj.diagramlist2] , 'String' , diagrams);
            
            sel1 = commodel.getSelectedDiagram(1);
            sel2 = commodel.getSelectedDiagram(2);
            
            set(obj.diagramlist1 , 'Value' , sel1);
            set(obj.diagramlist2 , 'Value' , sel2);
             
        end
        
        function syncCurveList(obj , commodel)
            curvind1 = commodel.getSelectedSolution(1);
            list1 = commodel.getSolutionList(1);
            buildList(obj.curvelist1 , list1 , curvind1);
            
            curvind2 = commodel.getSelectedSolution(2);
            list2 = commodel.getSolutionList(2);
            buildList(obj.curvelist2 , list2 , curvind2);            
        end
    end
    
end

function diagramSelectCallback(handle , commodel)
    ind = get(handle , 'UserData');
    valind = get(handle , 'Value');
    commodel.selectDiagram(ind, valind);
end

function curveSelectCallback(handle,commodel) 
    ind = get(handle , 'UserData');
    valind = get(handle , 'Value');
    commodel.selectSolution(ind, valind);
end

function setListContextMenu(handle , commodel , index)
    menu = uicontextmenu;
    GUIMenuItem(menu , @() commodel.viewCurve(index) , 'View' , commodel , @()  commodel.isValid(index)  , 'diagramSelectionChanged');

    plotm = uimenu(menu , 'Label' , 'Plot');
    
    GUIMenuItem(plotm , @() commodel.startSolutionPlot(index , 2)  , '2D', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged');    
    GUIMenuItem(plotm , @() commodel.startSolutionPlot(index , 3)  , '3D', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged');    

    
    
    GUIMenuItem(menu , @() commodel.renameSolution(index) , 'Rename', commodel , @()  commodel.isValid(index)  , 'diagramSelectionChanged' , 'Separator' , 'on');
    GUIMenuItem(menu , @() commodel.deleteSolution(index) , 'Delete', commodel , @()  commodel.isValid(index)  , 'diagramSelectionChanged');
    set(handle, 'UIContextMenu' , menu);

end

function setDiagramContextMenu(handle , commodel , index)
    menu = uicontextmenu;
   
    
    plotm = uimenu(menu , 'Label' , 'Plot');    
    GUIMenuItem(plotm , @() commodel.startDiagramPlot(index , 2)  , '2D', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged');
    GUIMenuItem(plotm , @() commodel.startDiagramPlot(index , 3)  , '3D', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged');    
    
    GUIMenuItem(menu , @() commodel.newDiagram(index)  , 'New', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged' , 'Separator' , 'on');
    GUIMenuItem(menu , @() commodel.renameDiagram(index)  , 'Rename', [] , @() commodel.getNrDiagrams() > 0  , 'diagramListChanged');
    GUIMenuItem(menu , @() commodel.deleteDiagram(index) , 'Delete', commodel , @() commodel.getNrDiagrams() > 1  , 'diagramListChanged');
    
    set(handle, 'UIContextMenu' , menu);
end

function diagramSelectionChanged(listhandle , commodel , event)
    ind = get(listhandle, 'UserData');
    
    if(event.index == ind)
        if (event.curveindex == event.ADJUSTINDEX)
            curveindex = get(listhandle, 'Value');
            len = length(event.curvelist);
            if (len == 0)
                %donothing
            elseif (curveindex > len)
               curveindex = len; 
               commodel.selectSolution(ind, curveindex);
            end
            buildList(listhandle , event.curvelist , curveindex);
            
        else
            buildList(listhandle , event.curvelist , event.curveindex);
        
        end
        
    end

end

function buildList(listhandle , list , listind)
    if ((listind > 0) && (~isempty(list)))
       set(listhandle , 'String' , list , 'Value' , listind, 'Enable' , 'on' , 'Selected' , 'off');
    else
        set(listhandle , 'String' , {' '} , 'Value' , 1 , 'Enable' , 'off');
    end

end
function keypress(handle , event, commodel, sideindex)
        if strcmp(event.Key, 'delete')
            commodel.deleteCurve(sideindex);
        end
end



