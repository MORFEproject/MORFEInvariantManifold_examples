classdef GUIMonitorPanel < handle
    
    properties
        panelhandle
        layoutstructure
        
        outputhandle
        statushandle
        
        session
    end
    
    methods
        function obj = GUIMonitorPanel(session, windowmanager, varargin)
            global sOutput
            obj.session = session;
            parent = windowmanager.demandWindow('monitor');
            pos = get(parent, 'Position');
            
            obj.panelhandle = uipanel(parent, 'Unit' , 'Pixels' ,  'BackgroundColor' , [0.85 0.85 0.85], 'DeleteFcn' , @(o,e) obj.destructor() ...
                , 'Position' , [0 0 pos(3) pos(4)] , varargin{:});
            
            
            mainnode = LayoutNode(-1 , -1 , 'vertical');
            
            mainnode.addNode( LayoutNode(1,1));
            obj.statushandle = uicontrol(obj.panelhandle,'Style' , 'edit' , 'String' , 'Ready' , 'BackgroundColor' , DefaultValues.BACKGROUNDCOLOR, 'Unit' , 'Pixels' ,...
                'HorizontalAlignment' , 'center' , 'Enable' , 'inactive', 'FontSize', 13 , 'FontWeight' , 'bold' );
            mainnode.addHandle(4,1,obj.statushandle , 'minsize',[Inf,Inf]);
            mainnode.addNode( LayoutNode(2,1));
            mainnode.addGUIobject(4,1, GUIPauseResumeButton( obj.panelhandle , @() obj.launchCurrentCurve()),'minsize',[Inf,Inf]);
            mainnode.addGUIobject(4,1, GUIStopCloseButton( obj.panelhandle , @() close(parent)),'minsize',[Inf,Inf]);
            mainnode.addNode( LayoutNode(1,1));
            
            
            %obj.outputhandle = uicontrol(obj.panelhandle ,'Style' , 'listbox', 'HorizontalAlignment' , 'left' , 'Enable' , 'inactive' , 'BackgroundColor' , 'white');
            obj.outputhandle = []; %disabled
            %mainnode.addHandle(40,1,obj.outputhandle , 'minsize',[Inf,Inf]);
            
            obj.layoutstructure = mainnode;
            mainnode.makeLayoutHappen( get(obj.panelhandle,'Position') );
            
            
            sOutput.setStatusHandle( obj.statushandle );
            sOutput.setOutputHandle( obj.outputhandle );
            set(obj.outputhandle , 'String' , {'Computational output:' , ''});
            
            LayoutNode.normalize(obj.panelhandle);
            
            contextmenu = uicontextmenu('Parent', parent);
            uimenu(contextmenu, 'Label' , 'Pause at special points (this continuation only)' , 'Callback', @(o,e) sOutput.setPauseSpecial());
            uimenu(contextmenu, 'Label' , 'Pause at each point (this continuation only)' , 'Callback', @(o,e) sOutput.setPauseAlways());
            uimenu(contextmenu, 'Label' , 'Disable pausing (this continuation only)' , 'Callback', @(o,e) sOutput.setPauseNever());
            set(obj.panelhandle, 'UIContextMenu', contextmenu);
            
            set(parent, 'Visible' , 'on', 'KeyPressFcn', @(o, e) pushkeydown(e,o,obj.panelhandle),'KeyReleaseFcn', @(o, e) releasekey(e,o,obj.panelhandle)); %, 'CloseRequestFcn' , @(o,e) obj.closefunction(session)); TODO FIXME TMP
            drawnow;
            
        end
        
        function destructor(obj)
            global sOutput
            sOutput.setStatusHandle([]);
            sOutput.setOutputHandle([]);
            
            obj.layoutstructure.destructor();
            delete(obj);
        end
        
        function launchCurrentCurve(obj)
            CLBrowsers.switchBrowser(obj.session);
        end
        
        
        function closefunction(obj , session)
            
            if (session.isUnlocked())
                %sOutput.performStop();
                session.getWindowManager().closeWindow('monitor');
            end
            
        end
        
        

        
        
    end
    
end

function pushkeydown(event, object, panel)
    if ~strcmp(event.Key, 'control')
        SessionOutputKeypress(object, event);
        return;
    end
    global sOutput;
    panel.UserData = panel.BackgroundColor;
    object.UserData = sOutput.pauseoption;
    sOutput.pauseoption = sOutput.NEVER;
    panel.BackgroundColor = [0.8 0.9 1];
    drawnow;
    sOutput.unpause();
end

function releasekey(event, object, panel)
    if ~strcmp(event.Key, 'control'); return; end
    global sOutput;
    
    panel.BackgroundColor = panel.UserData;
    sOutput.pauseoption = object.UserData;
    sOutput.pause()
end

