classdef GUIDiagramOrganizerModel < handle
    
    properties
        diagramlist;
        
        session
        sessioncurrentdiagramname
        
        diagramselection
        
        solutionselection
        
        system
        
        diagram = {[] , []};
    end
    
    events
        diagramSelectionChanged
        diagramListChanged
        curveSelection
    end
    
    methods
        function obj = GUIDiagramOrganizerModel(session , system)
            obj.session = session;
            obj.sessioncurrentdiagramname = session.getDiagram().getName(); %there is always a diagram
            
            obj.system = system;
            obj.reconfigure();

            
            if (length(obj.diagramlist) > 1)
                obj.diagramselection = [1 2];
            else
                obj.diagramselection = [1 1];
            end
            
            
            obj.solutionselection = [1 1];
 
        end
        

        function reconfigure(obj)
           obj.diagramlist = CLDiagram.getDiagramList(obj.system.getDiagramPath());
           obj.notify('diagramListChanged');
        end
        
        
        
        function selectDiagram(obj, index , diagramindex)
            obj.diagramselection(index) = diagramindex;
            obj.diagram{index} = [];
            obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(index, obj.getSolutionList(index)  , 1 ));
        end
        
        function selectSolution(obj, index , solutionindex)
            obj.solutionselection(index) = solutionindex;
            obj.notify('curveSelection' , GUIDiagramOrganizerModelEvent(index, [] , solutionindex ));
        end
        
        function solutionname = getSelectedSolution(obj, index)
            solutionname = obj.solutionselection(index);
        end
        
        function curvedescr = getCurve(obj , index)
            i = obj.solutionselection(index);
            cm = obj.getDiagram(index);
            curvedescr = cm.getCurveByIndex(i);
        end
        
        function bool = isValidSolution(obj, index)
            i = obj.solutionselection(index);
            bool =    (i > 0) && ( i <= length( obj.getSolutionList(index) ));
        end
        
        
        
        function list = getSolutionList(obj, index)
            list = CLDiagram.getSolutionList(fullfile(obj.system.getDiagramPath(), obj.diagramlist{obj.diagramselection(index)}));
        end
        
        
        function dindex = getSelectedDiagram(obj , index)
            dindex = obj.diagramselection(index);
        end
        
        function list = getDiagramList(obj)
            list = obj.diagramlist;
        end
        
        function cm = getDiagram(obj , index)
            if (isempty(obj.diagram{index}))
                obj.diagram{index} = CLDiagram(obj.system.getDiagramPath(), obj.diagramlist{obj.diagramselection(index)});
            end
            cm =  obj.diagram{index};
        end
        
        
        function moveCurve(obj, fromIndex, toIndex)
            fromCM = getDiagram(obj , fromIndex);
            toCM   = getDiagram(obj, toIndex);
            fromList = getSolutionList(obj, fromIndex);
            cname = fromList{obj.solutionselection(fromIndex)};
            

            
            toNames = toCM.getSolutionNames();
            newname = cname;
            while any(strcmp(newname, toNames))
               newname = inputdlg('Enter new name', 'name already exists in other diagram', 1, {newname});
               if isempty(newname); return; end
               newname = newname{1};
            end
            
            [status, message] = copyfile(fullfile(fromCM.getPath(), cname), fullfile(toCM.getPath(), newname));
             if (~status)
               errordlg(message, 'move curve'); 
             end           
            
            fromCM.deleteSolution(obj.session, cname(1:end-4)); %WARNING FIXME TODO with .mat
            
            obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(1, obj.getSolutionList(1)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));
            obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(2, obj.getSolutionList(2)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));
        end
        
        function moveFromLeft(obj)
           obj.moveCurve(1 , 2); 
            
        end
        function moveFromRight(obj)
           obj.moveCurve(2 , 1); 
        end
        
        function b = moveFromLeftValid(obj)
            b = (obj.diagramselection(1) ~= obj.diagramselection(2)) && ~isempty(obj.getSolutionList(1)); 
        end
        
        function b = moveFromRightValid(obj)
            b = (obj.diagramselection(1) ~= obj.diagramselection(2)) && ~isempty(obj.getSolutionList(2));             
        end
        function b = isValid(obj, index)
            b =  ~isempty(obj.getSolutionList(index));
        end
        
        
        function deleteCurve(obj, index)
            cmanager = getDiagram(obj , index);
            list = getSolutionList(obj, index);
            cname = list{obj.solutionselection(index)};

            if strcmp( questdlg(['Are you sure you want to delete "'   obj.diagramlist{obj.diagramselection(index)} '/'   cname  '"'],'Delete curve'), 'Yes')
                cmanager.deleteSolution(obj.session, cname(1:end-4));
                obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(1, obj.getSolutionList(1)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));
                obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(2, obj.getSolutionList(2)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));

            end
            
        end
        function renameSolution(obj , index)
            
            cmanager = getDiagram(obj , index);
            list = getSolutionList(obj, index);
            cname = list{obj.solutionselection(index)};
            newname = inputdlg('Enter new name' , 'Rename' , 1 , {cname});
            if (~isempty(newname))
                newname = newname{1};
                [s, mess] = cmanager.renameSolution(obj.session, cname(1:end-4) , newname(1:end-4)); %FIXME MAT FILE.
                if (~s)
                    errordlg(mess, 'error');
                end
                obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(1, obj.getSolutionList(1)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));
                obj.notify('diagramSelectionChanged' , GUIDiagramOrganizerModelEvent(2, obj.getSolutionList(2)  , GUIDiagramOrganizerModelEvent.ADJUSTINDEX  ));
            end
        end

        function renameDiagram(obj,index)
            oldname = obj.diagramlist{obj.diagramselection(index)} ;
            newname = inputdlg('Enter new diagramname' , 'diagramname' , 1 , {oldname});
            if (~isempty(newname))
            if (ismember(newname{1} , obj.diagramlist))
                errordlg('diagramname already exists' , 'Error');
            else
                systempath = obj.system.getDiagramPath();
                
                [s,mess] = movefile( fullfile(systempath, oldname) , fullfile(systempath, newname{1}) ,'f'); %FILES
                if (~s)
                    errordlg(mess, 'error');
                else
                    cm = obj.session.getDiagram();
                    if (~isempty(cm) &&  strcmp(cm.getPath() , fullfile(systempath, oldname)))
                        obj.session.solutionhandle.diagram.path = fullfile(systempath, newname{1});
                        obj.session.solutionhandle.diagram.name = newname{1};
                        obj.session.notify('solutionChanged');
                    end
                    
                end
                obj.reconfigure();
            end
            end
        end
        
        function deleteDiagram(obj,index)
            dname = obj.diagramlist{obj.diagramselection(index)} ;
            systempath = obj.system.getDiagramPath();
            if strcmp( questdlg(['Are you sure you want to delete the directory "'  dname  '"'],'Delete diagram'), 'Yes')
                [s, mess] = rmdir(fullfile(systempath, dname),'s');  %FILES
                
                if (~s)
                    errordlg(mess, 'error');
                else
                    cm = obj.session.getDiagram();
                    
                    
                    if (~isempty(cm) &&  strcmp(cm.getPath() , fullfile(systempath, dname)))
                        obj.session.changeSystem(obj.system); %force reloading of system.
                    end
                    
                    if (obj.diagramselection(1) > (length(obj.diagramlist) - 1)) %lists not yet updated
                        obj.diagramselection(1) = length(obj.diagramlist) - 1;
                    end
                    if (obj.diagramselection(2) > (length(obj.diagramlist) - 1)) %lists not yet updated
                        obj.diagramselection(2) = length(obj.diagramlist) - 1;
                    end                    
                    
                end
                obj.reconfigure();
                
            end
        end
        
        function newDiagram(obj, index)
            
           answer = inputdlg('Enter new diagramname' , 'diagramname');
           if (~isempty(answer))
                if (ismember(answer{1} , obj.diagramlist))
                    errordlg('diagramname already exists' , 'error');
                else
                    mkdir(fullfile(obj.system.getDiagramPath(), answer{1}));
                    obj.reconfigure();
                end
           end            
        end
        function nr = getNrDiagrams(obj)
           nr = length(obj.diagramlist); 
        end
      
        function startDiagramPlot(obj, index , dimension)
            list = obj.getSolutionList(index);
            listindex = obj.solutionselection(index);
            if listindex > 0 && listindex <= length(list)
                name = list{obj.solutionselection(index)};
                diagr = obj.getDiagram(index);
                CLDiagramPlotter(diagr, dimension, name);
            end
            
        end
        
        function [solution, name] = loadSelectedSolution(obj, index)
            list = obj.getSolutionList(index);
            listindex = obj.solutionselection(index);
            if listindex > 0 && listindex <= length(list)
                name = list{obj.solutionselection(index)};
                solutionpath = fullfile(obj.getDiagram(index).getPath(), name);
                solution = CompSolution.load(solutionpath);
            else
               name = ''; solution = []; 
            end
                
        end
        
        
        function startSolutionPlot(obj , sideindex , dimension)
            solution = obj.loadSelectedSolution(sideindex);
            if isempty(solution); return; end
            
            if dimension == 2
                plot(solution);
            else
                plot3(solution);
            end
        end
        
        
        function viewSolution(obj, sideindex)
            [solution, name] = obj.loadSelectedSolution(sideindex);
            if isempty(solution); return; end 
            solution.showTable(name);
            
        end
        
        
        
    end

    
    methods(Static)
        function cm = manageCurrentDiagram(session)
            sys = session.getSystem();
            cm = GUIDiagramOrganizerModel(session , sys);
            session.windowmanager.closeWindow('browser');
            fig = session.windowmanager.demandWindow('diagramorganizer');
            GUIDiagramOrganizerPanel(fig, cm);
            set(fig, 'Visible' , 'on', 'WindowStyle' , 'modal');
            
            %delete(fig);%FIXME REMOVE
        end
        
    end
end


