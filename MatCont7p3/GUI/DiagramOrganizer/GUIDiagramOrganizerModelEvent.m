classdef GUIDiagramOrganizerModelEvent < event.EventData

    
    properties
        index;
        curvelist;
        curveindex;
    end
    properties(Constant)
       ADJUSTINDEX = -4; 
    end
    methods
        function obj = GUIDiagramOrganizerModelEvent(index , list , curveindex)
           obj.index = index; 
           obj.curvelist = list;
           obj.curveindex = curveindex;
        end
    end
    
end

