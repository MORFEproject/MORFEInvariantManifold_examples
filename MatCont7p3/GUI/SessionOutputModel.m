classdef SessionOutputModel < handle
    
    properties
       listenerlist = {}; 
    end
    
    events
        stateChanged
        settingChanged
    end
    
    methods
        function addToListeners(obj, le)
           obj.listenerlist{end+1} = le; 
        end
        
        
        function destructor(obj)
            for i = 1:length(obj.listenerlist)
                delete (obj.listenerlist{i});
            end
            delete(obj);
        end
    end
    

end

