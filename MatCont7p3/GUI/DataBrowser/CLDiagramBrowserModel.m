classdef CLDiagramBrowserModel < handle
    
    
    properties
        index
        diagrampath
        list
        label
        
        passon;
        previewdata = [];
    end
    
    methods
        
        function obj = CLDiagramBrowserModel(passon , path,label)
            obj.label = label;
            obj.diagrampath = path;
            
            obj.passon = passon;
            
            obj.initialise();
        end
        
        function keypress(obj,e)
            if strcmp(e.Key, 'delete')
                obj.delete_call()
            end
        end
        function initialise(obj)
            obj.list = CLDiagram.getDiagramList(obj.diagrampath);
            if (obj.passon.preview)
                obj.previewdata = cell(1,length(obj.list));
                for i = 1:length(obj.list)
                    obj.previewdata{i} = CLDiagram(obj.diagrampath, obj.list{i});
                end
            end

            obj.initIndex();
        end
        
        
        
        
        function initIndex(obj)
            if isempty(obj.list)
                obj.index = -1;
            else
                obj.index = 1;
            end            
            
            if ~isempty(obj.passon.session.getSystem())
                
                sh = obj.passon.session.getSolutionHandle();
  

                if strcmp(obj.passon.session.getSystem().getName(), obj.passon.system.getName()) && ~isempty(sh) && ~isempty(sh.diagram)
                    currentname = sh.diagram.getName();
                    if ~isempty(currentname)
                        obj.index = find(strcmp(obj.list, currentname));
                    end
                end
            end
            
            if isempty(obj.index) || obj.index < 1;  obj.index = 1;  end
            
            
        end
        
        function list = getList(obj)
            list = obj.list;
        end
        
        function bool = goUp(obj)
            bool = 1;
        end
        function b = isValidSelection(obj)
            b = (obj.index > 0);
        end
        
        function newobj = selectItem(obj,index)
            passon = obj.passon;
            passon.diagram = obj.list{index};
            newobj = CLSolutionBrowserModel(passon, fullfile(obj.diagrampath, obj.list{index}),  obj.list{index}, obj.previewdata{index});
        end
        
        
        function label = getLabel(obj)
            label = obj.label;
        end
        function flist = getOperations(obj)
            flist{1} = cmdStruct(@() 'New' , @() obj.create_new());
        end
        function flist = getItemOperations(obj)
            flist{1} = cmdStruct(@() 'Load'  ,    @()  obj.load_call());
            flist{2} = cmdStruct(@() 'Rename',    @()  obj.rename_call());
            flist{3} = cmdStruct(@() 'Delete',    @()  obj.delete_call());
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
        function type = getType(obj)
            type = 2;
        end
        function backOnTop(obj)
            obj.initialise();
        end
        
        
        
        
        function create_new(obj)
            answer = inputdlg('Enter new diagramname' , 'diagramname');
            if (~isempty(answer))
                if (ismember(answer{1} , obj.list))
                    errordlg('diagramname already exists' , 'error');
                else
                    mkdir(fullfile(obj.diagrampath,answer{1}));
                    obj.initialise();
                    obj.passon.listreloaded();
                end
            end
            
        end
        
        function rename_call(obj)
            if (obj.index > 0)
                oldname = obj.currentItemName();
                newname = inputdlg('Enter new diagramname' , 'diagramname' , 1 , {oldname});
                if (~isempty(newname))
                    if (ismember(newname{1} , obj.list))
                        errordlg('diagramname already exists' , 'error');
                    else
                        [s, mess, ~] = movefile(fullfile(obj.diagrampath, oldname), fullfile(obj.diagrampath, newname{1}) ,'f');
                        if (~s)
                            errordlg(mess, 'error');
                        else
                            session = obj.passon.session;
                            if session.hasDiagram() && strcmp(session.getDiagram().getDiagramsPath(),obj.diagrampath)  && strcmp(session.getDiagram().getName(), oldname)
                                session.solutionhandle.diagram = CLDiagram(obj.diagrampath, newname{1});
                                session.notify('solutionChanged');
                            end
                        end
                        obj.initialise();
                        obj.passon.listreloaded();
                    end
                end
            end
        end
        
        function load_call(obj)
            if (obj.index > 0)
                obj.passon.session.changeState(obj.passon.system, CLDiagram(obj.diagrampath, obj.currentItemName()), [], []);
                obj.passon.killswitch();
            end
        end
        
        function delete_call(obj)
            if (obj.index > 0)
                name =  obj.currentItemName();
                
                if strcmp( questdlg(['Are you sure you want to delete the directory "'  name  '"'],'Delete diagram'), 'Yes')
                    [s, mess, ~] = rmdir(fullfile(obj.diagrampath, name),'s');
                    if (~s)
                        errordlg(mess, 'error');
                    else
                        session = obj.passon.session;
                        if session.hasDiagram() && strcmp(session.getDiagram().getDiagramsPath(),obj.diagrampath)  && strcmp(session.getDiagram().getName(), name)
                            session.solutionhandle = CLSolutionHandler.fromSystem(obj.passon.system);
                            session.notify('solutionChanged');
                        end
                        
                        obj.initialise();
                        obj.passon.listreloaded();
                    end
                end
            end
        end
        
        function str = getToolTipString(obj)
            str = 'double-click to view computations';
        end
        
        function  pr = getPreviewRenderer(obj)
            pr = GUIPreviewRendererDiagrams();
        end
        
    end
    
end

function s = cmdStruct(label, cmd)
s = struct('label' , label , 'cmd' , cmd);
end




