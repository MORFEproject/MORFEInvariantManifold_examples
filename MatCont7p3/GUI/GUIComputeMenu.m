classdef GUIComputeMenu < handle
    
    properties
        handle
        eventlistener1
        eventlistener2
        eventlistener3
        eventlistener4
        menuitems
        
        isLocked = @() 0;
    end
    
    methods
        function obj = GUIComputeMenu(parent, type,  session, varargin)
            if(strcmp(type,'context'))
                obj.handle = uicontextmenu('Parent', parent, 'DeleteFcn', @(o, e) obj.destructor());
            else
                obj.handle = uimenu(parent, 'Label' , 'Compute'  , 'Enable' , 'off', 'DeleteFcn',@(o, e)  obj.destructor(), varargin{:});
            end
            
            obj.eventlistener1 = session.addlistener('computationChanged' , @(o,e) obj.configure(session));
            obj.eventlistener2 = session.addlistener('settingChanged'     , @(o,e) obj.syncValidity(session));
            obj.eventlistener3 = session.addlistener('solutionChanged'    , @(o,e) obj.syncValidity(session));
            obj.eventlistener4 = session.addlistener('lockChanged'    , @(o,e) obj.syncValidity(session));
            obj.configure(session);
            
            
            obj.isLocked = @() session.isLocked();
        end
        function destructor(obj)
            delete(obj.eventlistener1);
            delete(obj.eventlistener2);
            delete(obj.eventlistener3);
            delete(obj.eventlistener4);
            delete(obj);
        end
        
        
        function configure(obj, session)
            obj.menuitems = {};
            hs = allchild(obj.handle);
            if (~isempty(hs)), delete(hs); end
            
            alist = session.getComputeActions();
            
            if isempty(alist) || obj.isLocked()
                set(obj.handle, 'Enable', 'off');
                
            else
                set(obj.handle, 'Enable', 'on');

                for k = 1:length(alist)
                   if ~isempty(alist(k).function)
                        separator = k-1 >= 1 && isempty(alist(k-1).function); 
                        obj.menuitems{end+1} = uimenu(obj.handle, 'Separator', CLbool2text(separator), 'Label', alist(k).label,...
                            'CallBack', @(o, e) session.compute(k), 'UserData', alist(k).valid);
                   end
                end
                
                
                obj.syncValidity(session);
            end
        end
        
        function syncValidity(obj, session)
            for k = 1:length(obj.menuitems)
                item = obj.menuitems{k};
                set(item, 'Enable', CLbool2text(feval(item.UserData, session, session.getSolution())));
            end
            set(obj.handle, 'Enable', CLbool2text(~obj.isLocked())); 
            
        end
        
        
    end
    
end

