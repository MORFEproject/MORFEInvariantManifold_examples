classdef CLSolutionBrowserModel < handle

    properties
        index
        label

        previewdata
        passon
        
        list
        diagram
    end
    
    methods
        
        function obj = CLSolutionBrowserModel(passon, path, label, diagram)
            obj.passon = passon;
            obj.label = label;
            obj.diagram = diagram;

            obj.initList();
            obj.initIndex();
        end
        function initList(obj)
            obj.previewdata = obj.diagram.getSolutions();
            obj.list = obj.diagram.getCleanSolutionNames();
        end
        
        function backOnTop(obj)
        end  
        
        
        function keypress(obj,e)  
            if strcmp(e.Key, 'delete')
                obj.delete_call();
            end            
        end      
        
        function list = getList(obj)
              list = obj.list;
        end 
        
        function initIndex(obj)
            if isempty(obj.list)
                obj.index = -1;
            else
                obj.index = 1;  
                if ~isempty(obj.passon.session.getSystem())
                    
                    sh = obj.passon.session.getSolutionHandle();
                    browserdiagram = obj.passon.diagram;
                    if ~ischar(browserdiagram); browserdiagram = browserdiagram.getName(); end
                    if strcmp(obj.passon.session.getSystem().getName(), obj.passon.system.getName()) && ~isempty(sh) && ~isempty(sh.solution) && strcmp(browserdiagram, sh.diagram.getName())
                        currentname = sh.getSolutionName();
                        if ~isempty(currentname)
                            obj.index = find(strcmp(obj.list, currentname));
                        end
                    end
                end
                if isempty(obj.index) || obj.index < 1;  obj.index = 1;  end
                
                
            end
        end          
        
        
        
        function bool = goUp(obj)
           bool = 1; 
        end
        function b = isValidSelection(obj)
           b = (obj.index > 0); 
        end
        
        function newobj = selectItem(obj,index) 
            passon = obj.passon;
            passon.diagram = obj.diagram;
            newobj = CLSwitchBrowserModel(passon, obj.previewdata{index}, obj.list{index});
        end
        
        function label = getLabel(obj)
            label = obj.label;
        end        
        
        function flist = getOperations(obj)
            flist{1} = cmdStruct(@() 'New'  ,    @() obj.createNew() );
            
        end
        function flist = getItemOperations(obj)

            flist{1} = cmdStruct(@() 'Load'  , @() obj.loadCurve());
            flist{2} = cmdStruct(@() 'View'  , @() obj.previewPanel());
            flist{3} = cmdStruct(@() 'Rename', @() obj.rename_call());
            flist{4} = cmdStruct(@() 'Delete', @() obj.delete_call());

        end
        function flist = getInheritedOperations(obj)
           flist = {}; 
        end
        
        function name = currentItemName(obj)
            if (obj.index < 0)
                i = 1;
            else
                i = obj.index; 
            end
            name = obj.list{i};
        end
        
        function data = getPreviewData(obj)
            data = obj.previewdata;
        end
        
        function data = getSelectedPreviewData(obj) %gevaarlijker
            if (obj.passon.preview && (obj.index >= 1))
                data = obj.previewdata{obj.index};
            else
                data = [];
            end
        end

        
        function createNew(obj)
            obj.passon.session.changeState(obj.passon.system, obj.diagram, [], []);
            %create empty settings.
            obj.passon.session.computation = CompConf();  %clear selected branch/curve/computation
            settings = CLSettings();
            settings.installSetting('system', obj.passon.system);
            obj.passon.session.changeSettings(settings);
            
            obj.passon.session.notify('settingsChanged');
            obj.passon.killswitch();
        end
        
        
        function previewPanel(obj)
            solution = obj.previewdata{obj.index};
            [b, msg] = solution.showTable(obj.list{obj.index}, obj.passon.session);
            
            if (~b)
                errordlg(msg, 'error');
                
            end
  
        end
        
        
        function loadCurve(obj) 
            solution = obj.previewdata{obj.index};
            solution.name = obj.list{obj.index};
            obj.passon.session.changeState(obj.passon.system, obj.diagram, solution, []);
            obj.passon.killswitch();
        end
        
     
        
        function rename_call(obj)
            if (obj.index > 0)
                oldname = obj.currentItemName();
                newname = inputdlg('Enter new curvename' , 'curvename' , 1 , {oldname});
                if (~isempty(newname))
                    
                    [s, mess] = obj.diagram.renameSolution(obj.passon.session, oldname , newname{1});
                    if (~s)
                       errordlg(mess, 'error'); 
                    end

                    obj.initList();
                    obj.passon.listreloaded();
                end
            end
        end
        
        
        function delete_call(obj)
            if (obj.index > 0)
                name =  obj.currentItemName();

                if strcmp( questdlg(['Are you sure you want to delete the curve "'  name  '"'],'Delete curve'), 'Yes')
                    obj.diagram.deleteSolution(obj.passon.session, name);
                    obj.initList();
                    obj.initIndex();
                    obj.passon.listreloaded();
                end
            end
            
        end
        
        function pr = getPreviewRenderer(obj)
            pr = GUIPreviewRendererSolutions();
        end
        function str = getToolTipString(obj)
            str = 'double-click to view points';
        end        
    end
    
    
    
end


function s = cmdStruct(label, cmd)
    s = struct('label' , label , 'cmd' , cmd);
end
