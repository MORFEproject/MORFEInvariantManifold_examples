classdef GUIMainPanel < handle

    properties
        panelhandle
        session
        mainnode
         
         
        labellist = { 'Class' , 'System', 'Curve' , 'Point Type' , 'Curve Type' , 'Derivatives' , 'Diagram' };
        
        labelstruct;
       	eventlistener
       
   
    
    end
    
    methods
        function obj = GUIMainPanel(parentfigure , session , varargin)
            if isempty(parentfigure)
               parentfigure = figure('Visible' , 'off', 'tag', 'matcont' , 'NumberTitle' , 'off' , 'Name' , 'MATCONTODE-DEMO', 'MenuBar' , 'none' );
            end
            
            obj.panelhandle = uipanel(parentfigure , 'Unit' , 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'BorderType' , 'none' ,  varargin{:});
            set(obj.panelhandle,'DeleteFcn' , @(o,e) obj.destructor(), 'ResizeFcn' , @(o,e) obj.onResize());
            
            
            obj.session = session;
            
            obj.mainnode = LayoutNode(-1,-1,'vertical');
            obj.mainnode.setOptions('sliderthreshold' , [100 100] , 'panel' , obj.panelhandle);
    
            
	    
            obj.addHeaderCompartment(parentfigure);
            obj.addSystemCompartment(parentfigure);
            obj.addCurveCompartment(parentfigure);
            
            
            obj.mainnode.makeLayoutHappen( get(obj.panelhandle, 'Position'));
            obj.eventlistener = [ obj.session.addlistener('settingsChanged', @(srv,ev) syncValues(obj, obj.session)), ...
                                  obj.session.addlistener('solutionChanged', @(srv,ev) syncValues(obj, obj.session))];
            obj.syncValues(obj.session);
            
            
           
           
           session.setInteractiveOutput(true);
           set( obj.panelhandle , 'Units','normalize', 'UserData', session); %NOTE: must be last!
           obj.session.setInteractiveOutput(true);
            
          oldclose = parentfigure.CloseRequestFcn;
          parentfigure.CloseRequestFcn = @(o, e) askClose(o, e, oldclose);
            
        end
        

        
        function destructor(obj)
            obj.session.setInteractiveOutput(false);
            obj.session.shutdownGUI();
            delete(obj.eventlistener)
            obj.mainnode.destructor();
      	    delete(obj);
        end
        
        function onResize(obj)
           set(obj.panelhandle,'Units' , 'Pixels');
           obj.mainnode.makeLayoutHappen( get(obj.panelhandle, 'Position'));
           set(obj.panelhandle,'Units' , 'normalize');
        end
        function addHeaderCompartment(obj, ~) %~=parentfigure
            panel = GUIPanel(obj.panelhandle , 0 ,  [5 1] , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                'String' , obj.labellist(1)) , 'halign' , 'l');
            obj.labelstruct.class = uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' ,  'String' , '' );
            subnode.addHandle(1,2, obj.labelstruct.class ,'halign' , 'l','minsize',[Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);
            obj.mainnode.addGUIobject(1,1, panel , 'minsize' , [Inf,Inf], 'margin' , [2,1]);
        end
        function addSystemCompartment(obj, parentfigure)
            panel = GUIPanel( obj.panelhandle , 8 , 5 , 'Title' , 'Current System', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                'String' , obj.labellist(2)) , 'halign' , 'l');
            obj.labelstruct.system = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' , 'String' , '' );
            subnode.addHandle(1,2, obj.labelstruct.system,'halign' , 'l','minsize',[Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);
            
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                'String' , obj.labellist(6)) , 'halign' , 'l');
            obj.labelstruct.derinfo = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' ,  'String' , '' );
            subnode.addHandle(1,2, obj.labelstruct.derinfo ,'halign' , 'l','minsize',[Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);
            
                    
                        
                
            context = uicontextmenu('Parent', parentfigure);
            
            global testcon;
            testcon = context;
            
 
            GUIMenuItem(context,@(o,e) CLBrowsers.systemBrowser(obj.session) , @() 'Load New System',...
                [] , @() obj.session.isUnlocked() , '');
            
            GUIMenuItem(context, @(o,e) obj.call_new_system(obj.session), @() 'Create New System',...
                [] , @() obj.session.isUnlocked() , '', 'Separator' , 'on');
            

            
            GUIMenuItem(context ,@(o,e) obj.call_edit_system(obj.session), @() 'Edit Current System',...
                obj.session , @() obj.session.hasSystem() && obj.session.isUnlocked() , 'settingsChanged' );
            

            
            GUIMenuItem(context , @(o,e) obj.call_userfunctions(obj.session) , @() 'Manage Userfunctions' , obj.session,...
                @() obj.session.hasSystem() && obj.session.isUnlocked() , 'settingsChanged');
            
            
            
            
            
            set([panel.handle, obj.labelstruct.system ,  obj.labelstruct.derinfo], 'UIContextMenu' , context);
            obj.mainnode.addGUIobject(2,1, panel,'minsize' , [Inf,Inf], 'margin' , [2,1]);
        end
        
        
  
        
        
        
        
        
        function addCurveCompartment(obj, parentfigure)

            panel = GUIPanel( obj.panelhandle , 8 , 5 , 'Title' , 'Current Curve', 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
  
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                 'String' , 'Name'  ),'halign' , 'l');
            obj.labelstruct.curvelabel = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' , 'String' , '');
            subnode.addHandle(1,2, obj.labelstruct.curvelabel,'halign' , 'l','minsize',[Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);

            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                'String' , 'Diagram'  ),'halign' , 'l');
            obj.labelstruct.diagram = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' ,  'String' , '');
            subnode.addHandle(1,2, obj.labelstruct.diagram,'halign' , 'l','minsize',[Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);            
            
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                'String' , 'Initial Point Type'  ),'halign' , 'l');
            obj.labelstruct.pointtype = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' ,  'String' , '');
            subnode.addHandle(1,2, obj.labelstruct.pointtype,'halign' , 'l','minsize' , [Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);            
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                 'String' , 'Curve Type'  ),'halign' , 'l');
            obj.labelstruct.curvetype = uicontrol( panel.handle, 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR , ...
                'HorizontalAlignment' , 'left' ,  'String','' );
            subnode.addHandle(1,2, obj.labelstruct.curvetype,'halign' , 'l','minsize' , [Inf,DefaultValues.LETTERDIMENSION(2)]);
            panel.mainnode.addNode(subnode);  
            
            subnode = LayoutNode(1,1);
            subnode.addHandle(1,1, ...
                uicontrol( panel.handle , 'Style' , 'text' , 'Units', 'Pixels' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'HorizontalAlignment' , 'left' , ...
                 'String' , 'Initializer'  ),'halign' , 'l');
            
             
            listbox = GUISelectBranchListbox(panel.handle , obj.session, @(x) x.toString(), 'Units', 'Pixels' ,   'HorizontalAlignment' , 'left',    'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR);
            obj.labelstruct.branchselector = listbox.handle;
            subnode.addHandle(1,2, obj.labelstruct.branchselector,'halign' , 'l','minsize' , [Inf,DefaultValues.LETTERDIMENSION(2)]);
            
            panel.mainnode.addNode(subnode);  
            panel.mainnode.addNode(subnode);
            
            installCurveContextMenu(obj,[panel.handle ,obj.labelstruct.curvelabel  ,obj.labelstruct.diagram  , obj.labelstruct.pointtype ,  obj.labelstruct.curvetype ], parentfigure);

            obj.mainnode.addGUIobject(3 , 1 , panel,'minsize' , [Inf,Inf] , 'margin' , [2,1] );
        end
        
        
        function installCurveContextMenu(obj, handlelist, parentfigure)
            context = uicontextmenu('Parent', parentfigure);

             GUIMenuItem(context, @(o,e) CLBrowsers.switchBrowser(obj.session), @() 'View Curve',...
                 obj.session , @() obj.session.hasSolution() && obj.session.isUnlocked() , 'solutionChanged');
             
             
             GUIMenuItem(context ,@(o,e) CLBrowsers.solutionBrowser(obj.session), @() 'View Diagram',...
                 obj.session , @() obj.session.hasSystem() && obj.session.isUnlocked() , 'settingsChanged' );

             
             GUIMenuItem(context , @(o,e) GUIMainPanel.renamecurve(obj.session) , @() 'Rename Curve' , obj.session,...
                  @() obj.session.hasSolution() && obj.session.isUnlocked() , 'solutionChanged'  , 'Separator' , 'on');
              
             GUIMenuItem(context , @(o,e) GUIMainPanel.newdiagram(obj.session) , @() 'New diagram' , obj.session,...
                  @() obj.session.hasSystem() && obj.session.isUnlocked() , 'settingsChanged' );             
             
             GUIWindowLaunchButton(context , obj.session.getWindowManager() , 'uimenu', ...
                 @(fhandle) GUIRenderSettings(fhandle, obj.session, [DefaultValues.SECTIONID.CONTINUER, DefaultValues.SECTIONID.INTEGRATOR]), 'continuer', 'Separator' , 'on');
             
             GUIWindowLaunchButton(context , obj.session.getWindowManager() , 'uimenu', ...
                 @(fhandle)  GUIRenderSettings(fhandle, obj.session, DefaultValues.SECTIONID.STARTER)  , 'starter');
            
            set(handlelist, 'UIContextMenu' , context);
        end        
        
        
        
        function syncValues(obj , session)
            set(obj.labelstruct.class      , 'String', 'ODE'); %class
            
            sh = session.getSolutionHandle();
            
            if session.hasSystem()
                set(obj.labelstruct.system     , 'String', session.getSystem().getName()     ); %system
                set(obj.labelstruct.derinfo    , 'String', session.getSystem().getDerInfo()  ); %derrs
            else
                set(obj.labelstruct.system     , 'String', '<no system>');
                set(obj.labelstruct.derinfo    , 'String', '-' ); %derrs
            end
            set(obj.labelstruct.curvelabel , 'String', sh.getSolutionName());
            set(obj.labelstruct.diagram    , 'String', sh.getDiagramName());

            pt = session.getPointType();
            lbl = pt.getLabel();
            if isempty(lbl)
               txt = lbl;
            else
                txt = sprintf('%s (%s)', pt.getName(), lbl);
            end
            set(obj.labelstruct.pointtype  , 'String', txt); %Point Type
            set(obj.labelstruct.curvetype  , 'String', sprintf('%s (%s)', session.getCompConf().getName() , session.getCompConf().getSolutionLabel() )); %Curve type
            
            
            %}
        end        
       
    end
    methods(Static)
        
        function  newdiagram(session)
            system = session.getSystem();
            path = system.getDiagramPath();
            
            list = CLDiagram.getDiagramList(path);
            
            answer = inputdlg('Enter new diagramname' , 'diagramname');
            if (~isempty(answer))
                if (ismember(answer{1} , list))
                    errordlg('diagramname already exists' , 'error');
                else
                    mkdir(fullfile(path,answer{1}));
                    session.changeState(system, CLDiagram(path , answer{1}), [], []);
                end
            end
            session.notify('diagramChanged');
        end
        
        
        function renamecurve(session)
            sh = session.getSolutionHandle();
            oldname = sh.getSolutionName();
            newname = inputdlg('Enter new curvename' , 'curvename' , 1 , { oldname });
            
            if (~isempty(newname))
                [s, mess] = sh.diagram.renameSolution(session, oldname, newname{1}, []);
                if (~s)
                    errordlg(mess, 'error');
                end
            end
        end
        

        

        function call_new_system(session)
            system = SysGUI.gui_loader(@SysGUI.new);
            if ~isempty(system)
                session.changeSystem(system);
            end
            
        end
        
        function call_edit_system(session)
            
            system = SysGUI.gui_loader(@() SysGUI.edit(session.getSystem().getName()));
            if ~isempty(system)
                session.changeSystem(system);
            end
            
        end
        function call_userfunctions(session)
            system = SysGUI.gui_loader(@() SysGUI.userfunctions(session.getSystem().getName()));
            if ~isempty(system)
                session.settings.installSetting('system', system);
                session.changeSettings(session.settings);
            end
            
        end
        
    end
end




function askClose(figurehandle, event, closefunction)
    if strcmp(questdlg('Are you sure you want to close MatCont?', 'Close MatCont', 'Yes', 'No', 'Yes'), 'Yes')
       closefunction(figurehandle, event);
    end
end
