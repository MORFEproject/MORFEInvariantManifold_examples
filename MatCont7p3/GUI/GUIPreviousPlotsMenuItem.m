classdef GUIPreviousPlotsMenuItem < handle
    
    properties
        handle
        
        eventlistener
        eventlistener2
        type;
        
    end
    
    methods
        function obj = GUIPreviousPlotsMenuItem(parent, type, session, varargin)
            
            
            if(strcmp(type,'context'))
                obj.handle = uicontextmenu;
            else
                obj.handle = uimenu(parent, 'Label' , 'Previous'  , 'Enable' , 'off' , varargin{:});
            end
            obj.type = type;
            
            
            obj.eventlistener = session.addlistener('solutionChanged' , @(o,e) obj.configure(session));
            obj.eventlistener2 = session.addlistener('computationChanged' , @(o,e) obj.configure(session));
            obj.configure(session);
        end
        
        
        function configure(obj, session)
            hs = allchild(obj.handle);
            if (~isempty(hs)), delete(hs); end
            
            pp = session.settings.getSetting('previousplots');
            list = pp.getConfigs();
            
            
            len = length(list);
            
            if (~strcmp(obj.type,'context'))
                set(obj.handle , 'Enable' , CLbool2text(len ~= 0));
            end
            
            for i = 1:length(list)
                uimenu(obj.handle, 'Label' , list{i}.toString() , 'Callback' ,  @(o,e) session.outputmanager.launch(list{i}));
            end
            
            
        end
        
    end
    
end

