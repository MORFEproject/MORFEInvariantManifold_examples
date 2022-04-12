classdef CLSystemBrowserModel < handle
    
    properties
        index
        systempath
        list
        label
        previewdata = [];  %DO mem cleanup after
        
        passon
    end
    
    methods
        
        function obj = CLSystemBrowserModel(passon, path, label)
            obj.label = label;
            obj.systempath = path;
            obj.passon = passon;
            obj.setUp();
                     
        end
        
        function setUp(obj)
            obj.list = CLSystemBrowserModel.getSystemList(obj.systempath);
            if (obj.passon.preview)
                obj.previewdata = cell(1,length(obj.list));
                for i = 1:length(obj.list)
                   obj.previewdata{i} = CLSystem(fullfile(obj.systempath, [obj.list{i} '.mat']));
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
                currentname = obj.passon.session.getSystem().getName();
                obj.index = find(strcmp(obj.list, currentname));
                
            end
            
            if isempty(obj.index) || obj.index < 1;  obj.index = 1;  end
        end
        
        
        function reSet(obj)
            len = length(obj.list);
            if (len == 0)
                obj.index = -1;
            elseif (obj.index > len)
                obj.index = len;
            end            
        end
        
        
        function list = getList(obj)
              list = obj.list;
        end
        
        function bool = goUp(obj)
            %disp('begining of the line');
           bool = 0; 
        end

        function newobj = selectItem(obj,index) 
            %fprintf(2, 'selecting %g \n', index);
            
            passon = obj.passon;
            system = obj.getSelectedPreviewData();
            passon.system = system;
            newobj = CLDiagramBrowserModel(passon, system.getDiagramPath(), obj.list{index});
            
        end
        
        function name = getItemName(obj,index)
            name = obj.list{index};
        end
        
        function name = currentItemName(obj)
            if (obj.index < 0)
                i = 1;
            else
                i = obj.index; 
            end
            name = obj.list{i};
        end
        function b = isValidSelection(obj)
           b = (obj.index > 0); 
        end
        
        function label = getLabel(obj)
            label = obj.label;
        end
        
        function flist = getOperations(obj)
            flist{1} = cmdStruct(@() 'New' , @() obj.new_system_call());
        end
        
        function flist = getItemOperations(obj)
            flist{1} = cmdStruct(@() 'Load'  ,    @()  obj.load_system_call());
            flist{2} = cmdStruct(@() 'Delete'  ,    @()  obj.delete_system_call());
            flist{3} = cmdStruct(@() 'Edit'  ,    @()  obj.edit_system_call());
            flist{4} = cmdStruct(@() 'Add/Edit Userfunctions'  ,    @()  obj.edit_uf_call());
        end
        
        function flist = getInheritedOperations(obj)
           flist = {}; 
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

        
        function backOnTop(obj)
           %fprintf(2, '%s: %s\n', mfilename, 'backontop');
        end
        
        function pr = getPreviewRenderer(~)
            pr = GUIPreviewRendererSystems();
        end
        
        
        function new_system_call(obj)
            system = SysGUI.gui_loader(@SysGUI.new);
            if ~isempty(system)
                obj.passon.killswitch();
                obj.passon.session.changeSystem(system);
            end
            
        end
        
        function load_system_call(obj)
            sys = obj.previewdata{obj.index};
            obj.passon.session.changeSystem(sys);
            
            
            [windowrestore, plotconfs] = obj.passon.session.loadFromFile('local');
            windowrestore(obj.passon.session);
            obj.passon.session.outputmanager.loadPlotConfs(plotconfs);
            obj.passon.session.notify('computationChanged');
            
            obj.passon.killswitch();
        end
        
        function edit_system_call(obj)
            fprintf(2, 'Altering a system might make older computed data uninterpretable.\n')
            
            system = obj.previewdata{obj.index};
            system = SysGUI.gui_loader(@() SysGUI.edit(system.getName()));
            
             if ~isempty(system)
                obj.passon.killswitch();
                obj.passon.session.changeSystem(system);
                %when not killed:
                %   obj.setUp();
                %   obj.passon.listreloaded();
             end     
        end
        
        function edit_uf_call(obj)
            system = obj.previewdata{obj.index};
            system = SysGUI.gui_loader(@() SysGUI.userfunctions(system.getName()));
             if ~isempty(system)
                obj.passon.killswitch();
                obj.passon.session.changeSystem(system);
                %when not killed:
                %   obj.setUp();
                %   obj.passon.listreloaded();
             end                
        end
        
        
        function delete_system_call(obj)
            session = obj.passon.session;
            system = obj.previewdata{obj.index};
            path = system.getDiagramPath();
            
            if strcmp( questdlg(['Are you sure you want to delete the system "'  system.getName() '"'],'Delete system'), 'Yes')
                [s, mess, ~] = rmdir(path ,'s');
                if (~s)
                    errordlg(mess, 'error');
                else
                    delete([path '.mat']);
                    delete([path '.m']);
                    
                    %conflict check:
                    if session.hasSystem() && strcmp(session.getSystem().getName(), system.getName())
                        session.changeSystem(CLSystem());
                    end
                    
                    obj.list(obj.index) = [];
                    obj.previewdata(obj.index) = [];
                    obj.reSet();
                    obj.passon.listreloaded();
                end
                
            end
        end
        
        
        function str = getToolTipString(obj)
            str = 'double-click to view diagrams';
        end
        
        function keypress(obj,e)
            if strcmp(e.Key, 'delete')
                obj.delete_system_call()
            end          
        end
    end
    
    methods(Static)
        function syslist = getSystemList(syspath)
            files = dir(syspath);
            allnames = {files.name};
            %filter files that end with .mat.
            names = allnames(cellfun(@(x) ~isempty(regexp(x, '\.mat$', 'ONCE')), allnames));
            %check if .m exists: cut of 'at' at the end of .mat'
            names = names(cellfun(@(x) ismember(x(1:end-2), allnames), names));
            syslist = cellfun(@(x) {x(1:end-4)}, names);  %cut of '.mat' at end

        end        
    end
    
end

function s = cmdStruct(label, cmd)
    s = struct('label' , label , 'cmd' , cmd);
end

