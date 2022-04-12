classdef GUISelectPointsMenu < handle
    
    
    properties
        handle
        
        eventlisteners = {};

    end
    
    methods
        function obj = GUISelectPointsMenu(parent, session, varargin)
            obj.handle = uimenu(parent, 'Label' , 'Initial Point' , varargin{:});
            set(obj.handle,'DeleteFcn' , @(o,e) obj.destructor());
            IP = CLSettingIP();
            
            
            n = length(IP.initiallist);
            for k = 1:n
                ip = IP.initiallist{k};
                if isobject(ip)
                    uimenu(obj.handle , 'Label' , ip.getName() ,  'Callback', @(o,e) session.changeInitPoint(ip.getLabel()) , varargin{:} , 'Separator' , CLbool2text( k-1>=1 && ~isobject(IP.initiallist{k-1})));
                end
            end
            delete(IP);
        end

        function destructor(obj)
            delete(obj);
        end
    end
    
end

