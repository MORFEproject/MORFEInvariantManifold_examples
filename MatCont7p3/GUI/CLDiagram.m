classdef CLDiagram
    
    properties
        path
        name
    end
    
    
    methods
        function obj = CLDiagram(diagramspath, diagramname)
            assert(exist(diagramspath, 'dir') ~= 0, 'diagramsdirectory does not exist')
            obj.path = fullfile(diagramspath, diagramname);
            
            if (exist(obj.path, 'dir') ~= 7)
                mkdir(obj.path);
            end
            
            obj.name = diagramname;
            
        end
        
        function p = getDiagramsPath(obj)
            [p, ~, ~] = fileparts(obj.path);
        end
        function p = getPath(obj)
           p = obj.path; 
        end
        
        function n = getName(obj)
           n = obj.name; 
        end
        
        function list = getSolutionNames(obj)
            list = obj.getSolutionList(obj.path);
        end
         function list = getCleanSolutionNames(obj)
            list = cellfun(@(x) {x(1:end-4)}, obj.getSolutionNames());  %remove .mat extension
        end       
        function path = suggestDefaultFile(obj, pointlabel, curvelabel, maxnumber)
            names = obj.getSolutionNames();
            pattern = [sprintf('^%s_%s', pointlabel, curvelabel), '\(\d+\).mat$'];
            
            k = length(names);
            while k >= 1 && isempty(regexp(names{k}, pattern, 'ONCE'))
               k = k - 1; 
            end
            if k < 1
                num = 1;
            else
                ii = find(names{k} == '(');
                num = str2double(names{k}(ii(end)+1:end-5)) + 1; %remove .mat
                if num > maxnumber
                   num = 1; 
                end
            end
            path = fullfile(obj.path, sprintf('%s_%s(%i).mat', pointlabel, curvelabel, num));
        end
        
        function solutions = getSolutions(obj)
            names = obj.getSolutionNames();
            paths = cellfun(@(x) fullfile(obj.path, x), names, 'UniformOutput', false);
            solutions = cellfun(@(x) CompSolution.load(x), paths, 'UniformOutput', false);
        end
        
        function [s, message] = renameSolution(obj, session, oldname, newname, solution)
            
            if ismember(newname, obj.getCleanSolutionNames())
                s = 0;
                message = 'Name already exists in Diagram';
                return;
            end
            [s, message] = movefile(fullfile(obj.path, [oldname '.mat']), fullfile(obj.path, [newname '.mat']) ,'f');
            
            if s
               session.pointloader.renaming(fullfile(obj.path, [oldname '.mat']), fullfile(obj.path, [newname '.mat']));
               
               if session.hasDiagram() && strcmp(session.getDiagram().getPath(), obj.getPath()) && strcmp(session.getSolutionHandle().getSolutionName(), oldname)
                   session.solutionhandle.solutionpath = fullfile(obj.path, [newname '.mat']);
                   session.notify('solutionChanged');
               end
               
            end
            
            

            
        end
        
        function deleteSolution(obj, session, name)
            
            %DELETE BUGFIX// this fixes some issues on windows systems.
            movefile(fullfile(obj.path, [name '.mat']), fullfile(obj.path, 'dummy') ,'f');
            delete(fullfile(obj.path, 'dummy')); 
            %%%% END BUGFIX
            
            if session.hasDiagram() && strcmp(session.getDiagram().getPath(), obj.getPath()) && strcmp(session.getSolutionHandle().getSolutionName(), name)
                session.changeState(session.getSystem(), obj, [], []);
            end
                
        end
        
        function labels = getSolutionLabels(obj, pointloader)
            names = obj.getSolutionNames();
            paths = cellfun(@(x) fullfile(obj.path, x), names, 'UniformOutput', false);
            if isempty(pointloader)
               labels =  cell(1, length(paths)); 
            else
               labels = cellfun(@(x) pointloader.getLabel(x), paths, 'UniformOutput', false);
            end
        end
        
        %{
        function plotdim(obj, dim)
        end
        function plot(obj);  obj.plotdim(2); end
        function plot3(obj);  obj.plotdim(3); end
        %}
    end

    methods(Static)
        function diagramlist = getDiagramList(path)
            files = sortfiles(dir(path));
            names = {files.name};
            names = names([files.isdir]);

            diagramlist = {};

            for i = 1:length(names)
                if ~strcmp(names{i}(1) , '.')
                    diagramlist{end+1} = names{i};
                end
            end
        end
        function list = getSolutionList(path)
            files = sortfiles(dir(path));
            list = {files.name};
            list = list(cellfun(@(x) ~isempty(regexp(x, '\.mat$', 'ONCE')), list));            
            
        end
        
    end
    
end

function files = sortfiles(files)
[~, ii] = sort([files.datenum]);
files = files(ii);
end

