classdef CLBrowsers
    methods(Static)
        

        
        function bm = systemBrowser(session)
            bm = CLBrowsers.constructBrowserModel(session, [], [], []);
            GUIBrowserPanel(session.windowmanager.demandWindow('browser'), session, bm, true);
            
        end

        function bm = solutionBrowser(session)
            bm = CLBrowsers.constructBrowserModel(session, session.getSystem(), session.getDiagram(), []);
            GUIBrowserPanel(session.windowmanager.demandWindow('browser'), session, bm, true);
            
        end
        function bm = diagramBrowser(session)
            bm = CLBrowsers.constructBrowserModel(session, session.getSystem(), [], []);
            GUIBrowserPanel(session.windowmanager.demandWindow('browser'), session, bm, true);
            
        end 
        
        function bm = switchBrowser(session)
           bm = CLBrowsers.constructBrowserModel(session, session.getSystem(), session.getDiagram(), session.getSolution());
           GUIBrowserPanel(session.windowmanager.demandWindow('browser'), session, bm, true);
        end
        
        function browsermodel = liteSwitchBrowser(session)  %Deprecated (FIXME TODO)
            if session.hasSolution()
                [browsermodel, passon] = constructModel(session);
                passon.system = session.getSystem();
                passon.diagram = session.getDiagram();
                browsermodel.initCurrentModel(CLSystemBrowserModel(passon, session.getSystemsPath(), 'Systems'))
                browsermodel.initCurrentModel(CLSwitchBrowserModel(passon, session.getSolution(),  session.getSolutionHandle().getSolutionName()))
                GUIBrowserPanel(session.windowmanager.demandWindow('browser'), session, browsermodel, false);
                %passon.solution = session.getSolution();
            end
            
        end
        
        
        function bm = constructBrowserModel(session, system, diagram, solution)
            [bm, passon] = constructModel(session);
            bm.current = CLSystemBrowserModel(passon, session.getSystemsPath(), 'Systems');
            %bm.current.index = 1;
            
            if ~isempty(system)
                passon.system = system;
                bm.stack{1} = bm.current;
                bm.current = CLDiagramBrowserModel(passon, system.getDiagramPath(), system.getName());
                %bm.current.index = 1;
                
                if ~isempty(diagram)
                    passon.diagram = diagram;
                    bm.stack{2} = bm.current;
                    bm.current = CLSolutionBrowserModel(passon, diagram.getPath(), diagram.getName(), diagram);
                    %bm.current.index = -1;
                    %if ~isempty(bm.current.list); bm.current.index = 1; end
                    
                    if ~isempty(solution)
                        passon.solution = solution;
                        bm.stack{3} = bm.current;
                        bm.current = CLSwitchBrowserModel(passon, solution, solution.name);
                        bm.current.index = 1;
                    end
                end
            end
        end
    end

end

function [browsermodel, passon] = constructModel(session)
browsermodel = CLBrowserModel();
passon = struct();
passon.session = session;
passon.preview = true;
passon.killswitch = @() browsermodel.shutDown() ;
passon.listreloaded = @() browsermodel.notify('listChanged');
passon.model = browsermodel;

end
