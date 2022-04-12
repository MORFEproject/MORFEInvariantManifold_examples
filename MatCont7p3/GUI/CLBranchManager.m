classdef CLBranchManager < handle
    properties(Constant, Hidden)
        COMPCONFIGS = CompConfigurations();
        
    end
    properties
        available;
        selected;
        
        prevpointlabel = '____';
        
    end
    events
        newlist
        selectionChanged
    end
    
    
    methods
        function obj = CLBranchManager(session)
            % FIXME FIMXE TODO FIXME TODO
            session.addlistener('initpointChanged', @(o, e) obj.configure(session)); %FIXME TODO destructor?
            obj.configure(session);
            
        end
        
        function configure(obj, session)
            pointlabel = session.settings.getValue('IP').getILabel();
            
            if ~strcmp(obj.prevpointlabel, pointlabel)
                obj.prevpointlabel = pointlabel;
                obj.available = obj.COMPCONFIGS.getConfigs(session.settings);
                obj.selected = 0;
                obj.notify('newlist');
                if ~isempty(obj.available)
                    obj.select(session, 1);
                end
            end
            
        end
        function sync(obj, session, computation, settings)
            pointlabel = settings.getValue('IP').getILabel();
            obj.prevpointlabel = pointlabel;
            obj.available = obj.COMPCONFIGS.getConfigs(settings);
            obj.selected = 0;
            
            for sel = 1:length(obj.available)
                if strcmp(obj.available{sel}.getLabel(), computation.getLabel())
                    obj.selected = sel;
                    obj.notify('newlist');
                    return;
                end
            end
            obj.selected = 0;
            if ~isempty(obj.available)
                fprintf(2, 'Solution loaded where computation does not match the initial point\n');
                obj.select(session, 1);
            end
            
        end
        
        
        
        function [l, index] = getCompConfigs(obj)
            l = obj.available;
            index = obj.selected;
        end
        
        function select(obj, session, index)
            obj.selected = index;
            cconf = obj.available{obj.selected};
            obj.notify('selectionChanged');
            session.selectCompConf(cconf);
        end
        
        function index = getSelectionIndex(obj)
            index = obj.selected;
        end
        function cconf = getSelection(obj)
           cconf = obj.available{obj.selected};
        end
        
        function b = validSelectionMade(obj)
            b = obj.selected > 0;
        end
        
        function setConf(obj, session, compconflabel)
            
            for k = 1:length(obj.available)
                if strcmp(compconflabel, obj.available{k}.getLabel())
                   obj.select(session, k);
                   return; 
                end
            end
        end
        
        
        % % % % % %
        
        function [labels, names, name] =  getConfVariants(obj, curvetype)
           labels = {};
           names = {};
           name = '';
           if strcmp(curvetype, 'O')
                confs = obj.COMPCONFIGS.OrbitConfs;
                labels = cellfun(@(x) x.getLabel(), confs, 'un', false);
                names = cellfun(@(x) x.getMethodName(), confs, 'un', false);
                name = 'Method';
           elseif strcmp(curvetype, 'ConnecHom')
                confs = obj.COMPCONFIGS.OrbitConfsHom;
                labels = cellfun(@(x) x.getLabel(), confs, 'un', false);
                names = cellfun(@(x) x.getMethodName(), confs, 'un', false);
                name = 'Method';               
           end
           %if strcmp(curvetype, 'ConnectHom/Het') FIXME TODO
           %end
           
        end
        
        function labels = getAllCurveLabels(obj)
            labels = obj.COMPCONFIGS.getCurveLabels();
        end
        
    end
    
end
