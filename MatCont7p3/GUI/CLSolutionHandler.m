classdef CLSolutionHandler < handle
    properties
       solution;
       solutionpath;
       diagram;
       
       
    end
    methods
        
        function n = getDiagramName(obj)
           if isempty(obj.diagram)
               n = '<no diagram selected>';
           else
               n = obj.diagram.getName();
           end
            
        end
        
        function n = getSolutionName(obj)
            if isempty(obj.solutionpath)
                if ~isempty(obj.solution)
                    n = '<nameless solution>';
                else
                    n = '<nothing>';
                end
            else
                [~, n, ~] = fileparts(obj.solutionpath); 
            end
        end
        
        function b = hasSolution(obj)
           b = ~isempty(obj.solution); 
        end
        
        function name = getDefaultName(obj, session)
           assert(~isempty(obj.diagram));
           settings = session.settings;
           name = obj.diagram.suggestDefaultFile(settings.getValue('IP').getLabel(), session.computation.getSolutionLabel(), settings.option_archive); 
        end
        
        function store(obj, path, solution)
           obj.solution = solution;
           obj.solutionpath = path;
           obj.solution.save(path);
            
        end
        
        
        function syncHandleWithMachine(obj)
            if ~isempty(obj.diagram) && ~isempty(obj.diagram.path)
                obj.diagram.path = correctSystemsDir(obj.diagram.path);
            end
            if ~isempty(obj.solutionpath)
                obj.solutionpath = correctSystemsDir(obj.solutionpath);
            end
            
        end
        
    end
    methods(Static)
        function obj = fromSystem(odesys)
                odesys.diagramInit();
                diagramspath = odesys.getDiagramPath();
                list = CLDiagram.getDiagramList(diagramspath);
                
            
                obj = CLSolutionHandler();
                if ~isempty(list)
                    diagramname = list{end};
                    obj.diagram = CLDiagram(diagramspath, diagramname);
                else
                    obj.diagram = [];
                end
                
                obj.solution = [];
                obj.solutionpath = [];
                
        end
        function obj = fromDiagram(diagram)
            obj = CLSolutionHandler();
            obj.diagram = diagram;
            obj.solution = [];
            obj.solutionpath = [];
        end
        
        function obj = fromSolution(diagram, solution)
            obj = CLSolutionHandler();
            obj.diagram = diagram;
            obj.solution = solution;
            obj.solutionpath = fullfile(diagram.path, [solution.name '.mat']);
        end
        
        function obj = fromSingleSolution(solution)
            obj = CLSolutionHandler();
            obj.solution = solution;
            
        end
        
        
    end
    
    
    

    
end
