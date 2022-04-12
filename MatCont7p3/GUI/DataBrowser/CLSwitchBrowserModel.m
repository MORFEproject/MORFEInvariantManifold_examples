classdef CLSwitchBrowserModel < handle
    
    
    properties
        index
        label
        passon
        
        solution
        switches
    end
    
    methods
        
        function obj = CLSwitchBrowserModel(passon, solution, label)
            obj.passon = passon;
            obj.label = label;
            obj.solution = solution;
            obj.switches = solution.listSwitches();
            obj.initIndex();
            
        end
        
        
        function backOnTop(obj)
        end
        
        function list = getList(obj)
            list = cellfun(@(k) {obj.switches.getName(k)}, num2cell(1:length(obj.switches)));
        end
        
        function bool = goUp(obj)
            bool = 1;
        end
        function b = isValidSelection(obj)
            b = (obj.index > 0);
        end
        
        function newobj = selectItem(obj, index)
            if index > 0
                newsettings = obj.switches.activateSwitch(index);
                if ~isempty(newsettings)
                    obj.passon.session.changeState(newsettings.system, obj.passon.diagram, obj.solution, newsettings);
                    obj.passon.killswitch();
                end
            end
            newobj = [];
        end
        function selectCurrentItem(obj)
           obj.selectItem(obj.index); 
        end
       
        
        function keypress(obj, e)
            %fprintf(2, 'keypress\n'); disp(e);
        end
        
        function label = getLabel(obj)
            label = obj.label;
        end
        
        
        function flist = getOperations(obj)
            flist{1} = cmdStruct(@() 'Load Curve' , @() obj.load_solution() );
            flist{2} = cmdStruct(@() 'View Settings' , @() GUIRenderSettings(obj.passon.session.getWindowManager().createFigure(), obj.solution.settings, [], false, ['Settings: ' obj.label]));  %[]: show everything (no filter), false: disallow editing.
            flist{3} = cmdStruct(@() 'View CurveData'  ,    @() obj.solution.showTable(obj.label, obj.passon.session));
            flist{4} = cmdStruct(@() 'Export'  ,    @() obj.export_solution());
            
        end
        
        function load_solution(obj)
            solution = obj.solution;
            solution.name = obj.label;
            obj.passon.session.changeState(obj.passon.system, obj.passon.diagram, solution, []);
            obj.passon.killswitch();
        end
        function export_solution(obj)
            varname = inputdlg('Export current solution to the MATLAB Command Window.', 'Enter variable name', 1, {'exported'});
            if isempty(varname); return; end
            varname = varname{1};
            
            assignin('base', varname, obj.solution);
            evalin('base', varname);
            
        end
        
        function initIndex(obj)
            if length(obj.switches) <= 0
                obj.index = -1;
            else
                obj.index = 1;
            end
        end
        function flist = getItemOperations(obj)
            flist{1} = cmdStruct(@() 'Select Point' ,    @()  obj.selectCurrentItem());
            
        end
        function flist = getInheritedOperations(obj)
            flist = {};
        end
        
        function data = getPreviewData(obj)
            data = obj.solution;
        end
        
        function data = getSelectedPreviewData(obj)
            if obj.index > 0
                data = obj.switches.getData(obj.index);
            else
               data = []; 
            end
        end

        function str = getToolTipString(obj)
            str = 'double-click to select a point';
        end
        function pr = getPreviewRenderer(obj)
            pr = obj.solution.getPreviewRenderer();
        end
    end
    
end

function s = cmdStruct(label, cmd)
s = struct('label' , label , 'cmd' , cmd);
end


